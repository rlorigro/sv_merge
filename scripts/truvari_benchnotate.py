import subprocess
import vcfpy
import sys
import os

from collections import OrderedDict


def benchnotate(truth_vcf_path, query_vcf_path, output_dir):
    args = [
        "truvari",
        "bench",
        "-b", truth_vcf_path,
        "-c", query_vcf_path,
        "--pctseq", str(0.9),
        "--sizemin", str(20),
        "-o", output_dir
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    fp_path = os.path.join(output_dir, "fp.vcf.gz")
    tp_path = os.path.join(output_dir, "tp-comp.vcf.gz")

    mapping = OrderedDict({
        'ID': 'TRUVARI_HG002_NIST',
        'Number': 1,
        'Type': 'Integer',
        'Description': 'Whether or not this variant was in the TP-comp.vcf of Truvari'
    })

    all_records = list()

    reader = vcfpy.Reader.from_path(fp_path)
    reader.header.add_info_line(mapping)

    for record in reader:
        reader.header.add_info_line(mapping)
        record.INFO["TRUVARI_HG002_NIST"] = 0
        all_records.append(record)

    reader = vcfpy.Reader.from_path(tp_path)
    reader.header.add_info_line(mapping)

    for record in reader:
        record.INFO["TRUVARI_HG002_NIST"] = 1
        all_records.append(record)

    out_path = os.path.join(output_dir, "input_annotated.vcf")

    chroms = {
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
    }

    writer = vcfpy.Writer.from_path(path=out_path, header=reader.header)

    for record in sorted(all_records, key=lambda x: (x.CHROM, x.POS)):
        if record.CHROM not in chroms:
            continue

        # print(record.INFO["TRUVARI_HG002_NIST"])

        writer.write_record(record)

    writer.close()

    args = [
        "bcftools",
        "view",
        "-I",
        out_path,
        "-O",
        "z",
        "-o",
        out_path + ".gz",
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None

    args = [
        "bcftools",
        "index",
        "-t",
        out_path + ".gz",
    ]

    sys.stderr.write(" ".join(args)+'\n')

    try:
        p1 = subprocess.run(args, check=True, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Status: FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
        sys.stderr.flush()
        return None


def main():
    truth_vcf_path = "/home/ryan/data/test_hapestry/vcf/hg002_nist/truvari_8x_bench/GRCh38_HG002-T2TQ100-V0.9_stvar.vcf.gz"
    vcf_path = "/home/ryan/data/test_hapestry/vcf/hg002_nist/truvari_8x_bench/HG002.truvari_collapsed.vcf.gz"
    output_dir = "/home/ryan/data/test_hapestry/vcf/hg002_nist/truvari_8x_bench/truvari_output"

    benchnotate(truth_vcf_path=truth_vcf_path, query_vcf_path=vcf_path, output_dir=output_dir)


if __name__ == "__main__":
    main()
