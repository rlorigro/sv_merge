import argparse
import math
import os


'''
name,length,is_ref,is_flank,coverage,identity,color
'''
def parse_nodes_csv(file_path):
    if not os.path.exists(file_path):
        return 0

    avg_coverage = 0

    column_name_to_index = dict()
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            data = line.strip().split(',')

            if i == 0:
                for index,item in enumerate(data):
                    column_name_to_index[item] = index

            else:
                nonref = data[column_name_to_index['is_ref']] == '0'

                if nonref:
                    x = data[column_name_to_index['coverage']]
                    avg_coverage += float(x) if "nan" not in x.lower() else 0

    if i > 0:
        avg_coverage /= i

    return avg_coverage


def parse_haps_csv(file_path):
    if not os.path.exists(file_path):
        return 0,0

    avg_coverage = 0
    avg_identity = 0

    column_name_to_index = dict()
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            data = line.strip().split(',')

            if i == 0:
                for index,item in enumerate(data):
                    column_name_to_index[item] = index

            else:
                a = data[column_name_to_index['coverage']]
                b = data[column_name_to_index['identity']]
                avg_coverage += float(a) if "nan" not in a.lower() else 0
                avg_identity += float(b) if "nan" not in b.lower() else 0

    if i > 0:
        avg_coverage /= i
        avg_identity /= i

    return avg_coverage, avg_identity


def parse_edges_csv(file_path):
    if not os.path.exists(file_path):
        return 0

    n_non_ref_edges_covered = 0
    n_non_ref_edges = 0

    column_name_to_index = dict()
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            data = line.strip().split(',')

            if i == 0:
                for index,item in enumerate(data):
                    column_name_to_index[item] = index

            elif i == 1:
                a = data[column_name_to_index['n_non_ref_edges_covered']]
                b = data[column_name_to_index['n_non_ref_edges']]

                n_non_ref_edges_covered += float(a) if "nan" not in a.lower() else 0
                n_non_ref_edges += float(b) if "nan" not in b.lower() else 0

            else:
                raise ValueError("edges.csv should only have 2 lines")

    avg_edges_covered = n_non_ref_edges_covered / (n_non_ref_edges if n_non_ref_edges > 0 else 1)

    return avg_edges_covered


'''
f(A,B) where A and B can be either a float or a tool name
'''
def compare(comparator_string, values_per_tool):
    tokens = comparator_string.strip().split(' ')

    if len(tokens) != 3:
        raise ValueError("Invalid comparator string has more than 3 tokens: " + comparator_string)

    a,comparator,b = tokens

    # if a is a tool name, get the value from the dict
    if a in values_per_tool:
        a = values_per_tool[a]
    else:
        try:
            a = float(a)
        except ValueError:
            raise Exception("Invalid comparator string f(A,B) has invalid A value: " + comparator_string + " with b = '" + a + "'")

    # if b is a tool name, get the value from the dict
    if b in values_per_tool:
        b = values_per_tool[b]
    else:
        try:
            b = float(b)
        except ValueError:
            raise Exception("Invalid comparator string f(A,B) has invalid B value: " + comparator_string + " with b = '" + b + "'")

    if comparator == ">":
        return a > b, a, b
    elif comparator == "<":
        return a < b, a, b
    elif comparator == ">=":
        return a >= b, a, b
    elif comparator == "<=":
        return a <= b, a, b
    elif comparator == "==":
        return a == b, a, b
    elif comparator == "!=":
        return a != b, a, b
    else:
        raise ValueError("Invalid comparator string f(A,B) has invalid comparator: " + comparator_string)


def write_results_to_bedgraph(output_path, results, name):
    with open(output_path, 'w') as out_file:
        # write header with data range from 0-1
        out_file.write("track type=bedGraph name=\"" + name + "\" autoScale=on\n")
        for dir_name,(a,b) in results.items():
            # parse dir name as BED format (TSV)
            tokens = dir_name.split('_')
            chrom = tokens[0]
            start,stop = list(map(int, tokens[1].split('-')))

            # compute diff on the values
            diff = abs(a - b)

            out_file.write(f"{chrom}\t{start}\t{stop}\t{diff}\n")

'''
Each directory has the following pattern:

[genomic_region]/
├── [tool_name]/
│    ├── alignments.gaf
│    ├── edges.csv
│    ├── graph.gfa
│    ├── haps.csv
│    ├── nodes.csv
│    ├── supported.vcf
│    └── unsupported.vcf

and there are multiple tools in the directory. The goal is to be able to filter the regions so that the ones matching 
the criteria can be written to a bed file. The criteria are provided by the user in the form of a comparator function,
e.g. --identity "hapestry > 0.9" or --identity "hapestry < truvari" 

There will be 4 possible stats that can be compared:
- hap_identity
- hap_coverage
- edges_covered
- nodes_covered

which can be obtained from the CSV files edges.csv nodes.csv and haps.csv

The format of edges.csv is:
n_alignments,n_edges,n_edges_covered,n_non_ref_edges,n_non_ref_edges_covered

The format of nodes.csv is:
name,length,is_ref,is_flank,coverage,identity,color

The format of haps.csv is:
name,length,is_ref,is_flank,coverage,identity

All CSVs contain header lines.

Nodes and haps CSVs need to be averaged within all rows of the file before comparison.
'''
def main(input_directory, output_directory, hap_identity, hap_coverage, edges_covered, nodes_covered):

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        raise ValueError("Output directory already exists: " + output_directory)

    hap_identity_results = dict()
    hap_coverage_results = dict()
    edges_covered_results = dict()
    nodes_covered_results = dict()

    # use os to iterate over the directories
    for d,dir_name in enumerate(os.listdir(input_directory)):
        if not os.path.isdir(os.path.join(input_directory, dir_name)):
            continue

        if d % 1000 == 0:
            print(d)

        hap_identity_per_tool = dict()
        hap_coverage_per_tool = dict()
        edges_covered_per_tool = dict()
        nodes_covered_per_tool = dict()

        try:
            # use os to iterate over the files in the directory
            for tool_name in os.listdir(os.path.join(input_directory, dir_name)):
                if not os.path.isdir(os.path.join(input_directory, dir_name, tool_name)):
                    continue

                # for each tool, extract the stats
                nodes_csv = os.path.join(input_directory, dir_name, tool_name, "nodes.csv")
                edges_csv = os.path.join(input_directory, dir_name, tool_name, "edges.csv")
                haps_csv = os.path.join(input_directory, dir_name, tool_name, "haps.csv")

                if hap_identity is not None or hap_coverage is not None:
                    avg_hap_coverage, avg_hap_identity = parse_haps_csv(haps_csv)
                    if hap_identity is not None:
                        hap_identity_per_tool[tool_name] = avg_hap_identity
                    if hap_coverage is not None:
                        hap_coverage_per_tool[tool_name] = avg_hap_coverage

                if edges_covered is not None:
                    avg_edges_covered = parse_edges_csv(edges_csv)
                    edges_covered_per_tool[tool_name] = avg_edges_covered

                if nodes_covered is not None:
                    avg_nodes_coverage = parse_nodes_csv(nodes_csv)
                    nodes_covered_per_tool[tool_name] = avg_nodes_coverage

            # use the comparator functions to see if this dir passes the criteria
            if hap_identity is not None:
                passing,a,b =  compare(hap_identity, hap_identity_per_tool)

                if passing:
                    hap_identity_results[dir_name] = (a,b)

            if hap_coverage is not None:
                passing,a,b =  compare(hap_coverage, hap_coverage_per_tool)

                if passing:
                    hap_coverage_results[dir_name] = (a,b)

            if edges_covered is not None:
                passing,a,b =  compare(edges_covered, edges_covered_per_tool)

                if passing:
                    edges_covered_results[dir_name] = (a,b)

            if nodes_covered is not None:
                passing,a,b =  compare(nodes_covered, nodes_covered_per_tool)

                if passing:
                    nodes_covered_results[dir_name] = (a,b)
        except Exception as e:
            print(e)
            print("ERROR: parsing failed for dir: " + dir_name)
            exit()


        # write the results to the bed file
        if len(hap_identity_results) > 0:
            name = "hap_identity" + " " + hap_identity
            write_results_to_bedgraph(os.path.join(output_directory, "hap_identity.bedgraph"), hap_identity_results, name)

        if len(hap_coverage_results) > 0:
            name = "hap_coverage" + " " + hap_coverage
            write_results_to_bedgraph(os.path.join(output_directory, "hap_coverage.bedgraph"), hap_coverage_results, name)

        if len(edges_covered_results) > 0:
            name = "edges_covered" + " " + edges_covered
            write_results_to_bedgraph(os.path.join(output_directory, "edges_covered.bedgraph"), edges_covered_results, name)

        if len(nodes_covered_results) > 0:
            name = "nodes_covered" + " " + nodes_covered
            write_results_to_bedgraph(os.path.join(output_directory, "nodes_covered.bedgraph"), nodes_covered_results, name)

    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare evaluation results')
    parser.add_argument('--input', help='Input directory')
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--identity', help='Identity comparison')
    parser.add_argument('--coverage', help='Coverage comparison')
    parser.add_argument('--edges_covered', help='Edges covered comparison')
    parser.add_argument('--nodes_covered', help='Nodes covered comparison')

    # input = "/home/ryan/data/test_hapestry/run/hprc_all_47_chr1/parameter_search_fix_case/full_eval/chr1_evaluation/"
    # output = "/home/ryan/data/test_hapestry/run/hprc_all_47_chr1/parameter_search_fix_case/full_eval/chr1_evaluation_compare"
    # hap_identity = "hapestry < truvari"
    # hap_coverage = None
    # edges_covered = None
    # nodes_covered = None

    main(
        parser.parse_args().input,
        parser.parse_args().output,
        parser.parse_args().identity,
        parser.parse_args().coverage,
        parser.parse_args().edges_covered,
        parser.parse_args().nodes_covered
    )

'''
EXAMPLE USAGE:
python3 compare_evaluation.py \
--input /home/ryan/data/test_hapestry/run/hprc_all_47_chr1/parameter_search_fix_case/full_eval/chr1_evaluation/ \
--output /home/ryan/data/test_hapestry/run/hprc_all_47_chr1/parameter_search_fix_case/full_eval/chr1_evaluation_compare \
--identity "hapestry < truvari" \
--coverage "hapestry < 0.9"
'''
