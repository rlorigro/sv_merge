from torch.utils.data.dataset import Dataset
import numpy as np
import torch

import vcfpy


def get_type_index(record: vcfpy.Record):
    if record.ALT[0].type == "INS":
        return 0
    elif record.ALT[0].type == "DEL":
        return 1
    elif record.ALT[0].type == "DUP":
        return 2
    else:
        return 3


def load_features_from_vcf(x: list,y: list, vcf_path: str, truth_info_name: str, annotation_name: str, filter_fn=None, contigs=None):
    print(vcf_path)
    reader = vcfpy.Reader.from_path(vcf_path)

    type_vector = [0,0,0,0]

    bp_norm = 15000

    for r,record in enumerate(reader):
        info = record.INFO
        ref_length = float(len(record.REF))/bp_norm
        alt_length = float(len(record.ALT[0].serialize()))/bp_norm

        if filter_fn is not None:
            if not(filter_fn(record)):
                continue

        if contigs is not None:
            if record.CHROM not in contigs:
                continue

        if record.calls is None:
            exit("ERROR: no calls in record: " + record.ID)
        elif len(record.calls) != 1:
            exit("ERROR: multiple calls in record: " + record.ID)

        # q = 0, p(correct) <= 0.0  Merged with above
        # q = 1, p(correct) <= 0.5  Merged with above
        # q = 2, p(correct) <= 0.75          i = 0 + 1
        # q = 3, p(correct) <= 0.875         i = 1 + 1
        # q = 4, p(correct) <= 0.9375        i = 2 + 1
        # q = 5, p(correct) <= 0.96875       i = 3 + 1
        # q = 6, p(correct) <= 0.984375      i = 4 + 1
        # q = 7, p(correct) <= 0.9921875     i = 5 + 1

        #          F NonSpan         R NonSpan           F Span           R Span         _ is_tandem
        #     |                 |                 |                 |                |  /  _ length of region evaluated
        # q    2  3  4  5  6  7  2  3  4  5  6  7  2  3  4  5  6  7  2  3  4  5  6  7  /  /
        # i 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
        # [ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, x]
        if truth_info_name.lower() == "hapestry":
            is_true = float(info["HAPESTRY_REF_MAX"]) > 0.9
        elif truth_info_name is not None:
            is_true = info[truth_info_name]
        else:
            is_true = 0

        caller_support = [0,0,0]
        caller_support[0] = info["SUPP_PAV"] if "SUPP_PAV" in info else 0
        caller_support[1] = info["SUPP_PBSV"] if "SUPP_PBSV" in info else 0
        caller_support[2] = info["SUPP_SNIFFLES"] if "SUPP_SNIFFLES" in info else 0

        t = type_vector
        t[get_type_index(record)] = 1

        y.append(is_true)
        x.append([])

        x[-1].extend(type_vector)
        x[-1].extend(caller_support)
        x[-1].append(float(ref_length)/float(bp_norm))
        x[-1].append(float(alt_length)/float(bp_norm))
        x[-1].append(record.QUAL if record.QUAL is not None else 0)
        # x[-1].append(float(info["SVLEN"][0])/float(bp_norm) if "SVLEN" in info else 0)
        x[-1].append(float(info["STDEV_POS"])/float(bp_norm) if "STDEV_POS" in info else 0)
        x[-1].append(float(info["STDEV_LEN"])/float(bp_norm) if "STDEV_LEN" in info else 0)

        reads = list(map(float,info["HAPESTRY_READS"]))
        reads[-1] /= bp_norm

        if annotation_name.lower() == "hapestry":
            x[-1].extend(reads)
            max_align_score = info["HAPESTRY_READS_MAX"]
            x[-1].append(max_align_score)

        elif annotation_name.lower() == "sniffles":
            call = record.calls[0]
            x[-1].append(call.data["GQ"])
            x[-1].append(call.data["DR"])
            x[-1].append(call.data["DV"])
            x[-1].extend(reads[-2:])


class VcfDataset(Dataset):
    def __init__(self, vcf_paths: list, truth_info_name, annotation_name, filter_fn=None, contigs=None):
        x = list()
        y = list()

        for p in vcf_paths:
            load_features_from_vcf(x=x,y=y,vcf_path=p,truth_info_name=truth_info_name,annotation_name=annotation_name,filter_fn=filter_fn,contigs=contigs)

        x = np.array(x)
        y = np.array(y)

        x_dtype = torch.FloatTensor
        y_dtype = torch.FloatTensor     # for MSE Loss

        self.length = x.shape[0]

        self.x_data = torch.from_numpy(x).type(x_dtype) + 1e-6
        self.y_data = torch.from_numpy(y).type(y_dtype)

    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]

    def __len__(self):
        return self.length

