from torch.utils.data.dataset import Dataset
import numpy as np
import torch
import sys

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


def load_features_from_vcf(x: list,y: list, feature_names: list, vcf_path: str, truth_info_name: str, annotation_name: str, filter_fn=None, contigs=None):
    print(vcf_path)
    reader = vcfpy.Reader.from_path(vcf_path)

    type_vector = [0,0,0,0]

    for r,record in enumerate(reader):
        info = record.INFO

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
            if truth_info_name not in info:
                sys.stderr.write("ERROR: truth info not found in record: " + str(record.ID) + "\n")
                continue

            is_true = info[truth_info_name]
        else:
            is_true = 0

        ref_length = float(len(record.REF))
        alt_length = float(len(record.ALT[0].serialize()))

        t = type_vector
        t[get_type_index(record)] = 1

        caller_support = [0,0,0]
        caller_support[0] = info["SUPP_PAV"] if "SUPP_PAV" in info else 0
        caller_support[1] = info["SUPP_PBSV"] if "SUPP_PBSV" in info else 0
        caller_support[2] = info["SUPP_SNIFFLES"] if "SUPP_SNIFFLES" in info else 0

        hapestry_data = list(map(float,info["HAPESTRY_READS"]))

        y.append(is_true)
        x.append([])

        q_score = record.QUAL if record.QUAL is not None else 0
        stdev_pos = float(info["STDEV_POS"]) if "STDEV_POS" in info else 0
        stdev_len = float(info["STDEV_LEN"]) if "STDEV_LEN" in info else 0

        x[-1].append(q_score)
        if r == 0:
            feature_names.append("q_score")

        x[-1].append(stdev_len)
        if r == 0:
            feature_names.append("stdev_len")

        x[-1].append(stdev_pos)
        if r == 0:
            feature_names.append("stdev_pos")

        x[-1].append(ref_length)
        if r == 0:
            feature_names.append("ref_length")

        x[-1].append(alt_length)
        if r == 0:
            feature_names.append("alt_length")

        x[-1].extend(type_vector)
        if r == 0:
            feature_names.extend(["type_" + str(i) for i in range(len(type_vector))])

        x[-1].extend(caller_support)
        if r == 0:
            feature_names.extend(["caller_support_" + str(i) for i in range(len(caller_support))])

        if annotation_name.lower() == "hapestry":
            max_align_score = info["HAPESTRY_READS_MAX"]

            x[-1].extend(hapestry_data)
            if r == 0:
                feature_names.extend(["hapestry_data_" + str(i) for i in range(len(hapestry_data))])

            x[-1].append(max_align_score)
            if r == 0:
                feature_names.append("max_align_score")

        elif annotation_name.lower() == "sniffles":
            call = record.calls[0]
            is_tandem = hapestry_data[-2]

            x[-1].append(call.data["GQ"])
            if r == 0:
                feature_names.append("GQ")

            x[-1].append(call.data["DR"])
            if r == 0:
                feature_names.append("DR")

            x[-1].append(call.data["DV"])
            if r == 0:
                feature_names.append("DV")

            x[-1].append(is_tandem)
            if r == 0:
                feature_names.append("is_tandem")

        elif annotation_name.lower() == "kanpig":
            is_tandem = hapestry_data[-2]
            call = record.calls[0]

            gt = list(map(float,call.data["GT"].split("|"))) if '|' in call.data["GT"] else [0,0]

            x[-1].extend(gt)
            if r == 0:
                feature_names.extend(["GT_" + str(i) for i in range(len(gt))])

            x[-1].extend(call.data["AD"])
            if r == 0:
                feature_names.extend(["AD_" + str(i) for i in range(len(call.data["AD"]))])

            x[-1].append(call.data["DP"])
            if r == 0:
                feature_names.append("DP")

            x[-1].append(call.data["GQ"])
            if r == 0:
                feature_names.append("GQ")

            x[-1].append(call.data["SQ"])
            if r == 0:
                feature_names.append("SQ")

            # x[-1].extend(call.data["ZS"])
            # if r == 0:
            #     feature_names.extend(["ZS_" + str(i) for i in range(len(call.data["ZS"]))])
            #
            # x[-1].extend(call.data["SS"])
            # if r == 0:
            #     feature_names.extend(["SS_" + str(i) for i in range(len(call.data["SS"]))])

            x[-1].append(is_tandem)
            if r == 0:
                feature_names.append("is_tandem")

        if r == 0:
            if len(feature_names) != len(x[-1]):
                print(feature_names)
                print("ERROR: feature names and data length mismatch: names:%d x:%d" % (len(feature_names), len(x[-1])))
                exit()


class VcfDataset(Dataset):
    def __init__(self, vcf_paths: list, truth_info_name, annotation_name, filter_fn=None, contigs=None):
        x = list()
        y = list()

        self.feature_indexes = None

        for p,path in enumerate(vcf_paths):
            feature_names = list()
            load_features_from_vcf(x=x,y=y,feature_names=feature_names,vcf_path=path,truth_info_name=truth_info_name,annotation_name=annotation_name,filter_fn=filter_fn,contigs=contigs)

            if p == 0:
                self.feature_indexes = {feature_names[i]: i for i in range(len(feature_names))}

        x = np.array(x)
        y = np.array(y)

        x_dtype = torch.FloatTensor
        y_dtype = torch.FloatTensor     # for MSE Loss

        self.length = x.shape[0]

        self.x_data = torch.from_numpy(x).type(x_dtype)
        self.y_data = torch.from_numpy(y).type(y_dtype)

        self.x_data = torch.nn.functional.normalize(self.x_data, dim=0)
        self.x_data += 1e-12

        self.filter_fn = filter_fn

    def __getitem__(self, index):
        return self.x_data[index], self.y_data[index]

    def __len__(self):
        return self.length

