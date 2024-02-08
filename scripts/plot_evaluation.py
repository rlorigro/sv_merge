from collections import defaultdict
import argparse
import numpy
import sys
import re
import os

from matplotlib import pyplot
import matplotlib

import numpy as np
from matplotlib.ticker import FormatStrFormatter


def plot_6_major_distributions(histograms: dict):

    axes_locations = {
        "alignment_identity_avg": (0,0),
        "haplotype_coverage_avg": (0,1),
        "n_alignments": (0,2),
        "nonref_nodes_fully_covered": (1,0),
        "nonref_edges_covered": (1,1),
        "nonref_nodes": (1,2)
    }

    is_reverse = {
        "nonref_alignment_identity_avg":True,
        "nonref_haplotype_coverage_avg":True,
        "cluster_alignment_identity_avg":True,
        "cluster_coverage_avg":True,
        "alignment_identity_avg":True,
        "haplotype_coverage_avg":True,
        "n_alignments":False,
        "nonref_nodes_fully_covered":True,
        "nonref_edges_covered":True,
        "nonref_nodes":False
    }

    # Show the full cumulative distribution for the following stats:
    figure,axes = pyplot.subplots(nrows=2, ncols=3)

    for s,(sample_label,series) in enumerate(histograms.items()):
        for series_label,coord in axes_locations.items():
            i,j = coord

            print(sample_label)
            h,bins = histograms[sample_label][series_label]
            print(h.shape)

            h = numpy.append(numpy.array([0]), h)
            bins = numpy.append(numpy.array([0]), bins)
            print(h.shape)

            # if is_reverse[series_label]:
            #     h = numpy.flip(h)

            cdf = numpy.cumsum(h)

            if is_reverse[series_label]:
                cdf = 1 - cdf

            print(bins[:10])
            print(bins[-10:])
            print(cdf[:10])
            print(cdf[-10:])
            print()
            axes[i][j].plot(bins,cdf, label=sample_label, alpha=0.7, linewidth=2)
            # axes[i][j].plot(bins,h, label=sample_label, alpha=0.7)
            axes[i][j].set_xlabel(series_label)

            axes[i][j].set_ylim([0,1.05])

            if is_reverse[series_label] and s == 0:
                axes[i][j].invert_xaxis()

    axes[0][0].legend()

    pyplot.show()
    pyplot.close()


def plot_all(histograms: dict):
    transposed_data = defaultdict(dict)

    for s,(sample_label,series) in enumerate(histograms.items()):
        for series_label,data in series.items():
            transposed_data[series_label][sample_label] = data

    is_reverse = {
        "nonref_alignment_identity_avg":True,
        "nonref_haplotype_coverage_avg":True,
        "cluster_alignment_identity_avg":True,
        "cluster_coverage_avg":True,
        "alignment_identity_avg":True,
        "haplotype_coverage_avg":True,
        "n_alignments":False,
        "nonref_nodes_fully_covered":True,
        "nonref_edges_covered":True,
        "nonref_nodes":False,
        "n_haplotype_clusters":False,
        "n_haplotypes":False,
        "n_nonref_haplotypes":False,
        "nonref_edges":False,
        "nonref_nodes_bps":False,
        "nonref_nodes_bps_covered":False,
        "nonref_nodes_partially_covered":False,
        "supported_del":True,
        "supported_dup":True,
        "supported_ins":True,
        "supported_inv":True,
        "supported_vcf_records":True,
        "sv_length_supported":False,
        "sv_length_unsupported":False,
        "unsupported_del":False,
        "unsupported_dup":False,
        "unsupported_ins":False,
        "unsupported_inv":False,
        "unsupported_vcf_records":False,
    }
    for s,[series_label,item] in enumerate(transposed_data.items()):
        if "n_alignments" not in series_label:
            continue

        fig = pyplot.figure()
        axes = pyplot.axes()

        all_series = list()

        for sample_label,series in item.items():
            print(sample_label)
            h,bins = series
            print(h.shape)

            h = numpy.append(numpy.array([0]), h)
            bins = numpy.append(numpy.array([0]), bins)
            print(h.shape)

            cdf = numpy.cumsum(h)

            if is_reverse[series_label]:
                cdf = 1 - cdf

            all_series.append(cdf)

            print(bins[:10])
            print(bins[-10:])
            print(cdf[:10])
            print(cdf[-10:])
            print()

            p = axes.plot(bins,cdf, label=sample_label, alpha=0.7, linewidth=0.8)
            c = p[-1].get_color()

            axes.plot(bins,h, alpha=0.7, linewidth=0.8, color=c)

            # axes.plot(bins,h, label=sample_label, alpha=0.7, linewidth=0.8)

            axes.set_xlabel(series_label)
            axes.set_ylim([0,1.05])

        if is_reverse[series_label]:
            axes.invert_xaxis()

        # all_series = numpy.stack(all_series)
        # print(all_series.shape)
        # x_variance = numpy.var(all_series, axis=0)
        # print(x_variance.shape)
        #
        # i_max_variance = numpy.argmax(x_variance)
        #
        # i_min = None
        # i_max = None
        #
        # for v,variance in enumerate(x_variance[i_max_variance:]):
        #     i_max = i_max_variance + v
        #     if variance < 0.5e-10:
        #         break
        #
        # for v,variance in enumerate(reversed(x_variance[:i_max_variance])):
        #     i_min = i_max_variance - v
        #     if variance < 0.5e-10:
        #         break
        #
        # if i_max > 0.8*len(bins):
        #     i_max = len(bins) - 1
        #
        # if i_min < 0.1*len(bins):
        #     i_min = 0
        #
        # x_max = bins[i_max]
        # x_min = bins[i_min]
        #
        # axes.axvline(x_max, color="gray")
        # axes.axvline(x_min, color="gray")
        #
        # axes.plot(bins,x_variance, label="Variance", alpha=0.7, linewidth=3)

        axes.legend()
        fig.set_size_inches(8,6)
        # fig.savefig(series_label + ".svg", dpi=200)
        pyplot.show()
        pyplot.close(fig)



def main(input_dir):
    # Arbitrary histogram bounds TODO: determine n_haps automatically
    n_haps = 95
    n_node_max = 300
    length_max = 8000

    data_ranges = {
        "alignment_identity_avg": [0,1],
        "cluster_alignment_identity_avg": [0,1],
        "cluster_coverage_avg": [0,1],
        "haplotype_coverage_avg": [0,1],
        "n_alignments": [0,n_node_max],
        "n_haplotype_clusters": [0,n_haps],
        "n_haplotypes": [0,n_haps],
        "n_nonref_haplotypes": [0,n_haps],
        "nonref_alignment_identity_avg": [0,1],
        "nonref_edges": [0,n_node_max],
        "nonref_edges_covered": [0,1],
        "nonref_haplotype_coverage_avg": [0,1],
        "nonref_nodes": [0,n_node_max],
        "nonref_nodes_bps": [0,length_max],
        "nonref_nodes_bps_covered": [0,1],
        "nonref_nodes_fully_covered": [0,1],
        "nonref_nodes_partially_covered": [0,1],
        "supported_del": [0,1],
        "supported_dup": [0,1],
        "supported_ins": [0,1],
        "supported_inv": [0,1],
        "supported_vcf_records": [0,1],
        "sv_length_supported": [0,length_max],
        "sv_length_unsupported": [0,length_max],
        "unsupported_del": [0,1],
        "unsupported_dup": [0,1],
        "unsupported_ins": [0,1],
        "unsupported_inv": [0,1],
        "unsupported_vcf_records": [0,1],
    }

    data = defaultdict(dict)


    # Get top level directories which correspond to each VCF that was evaluated
    for d,sample_label in enumerate(os.listdir(input_dir)):
        subdirectory_path = os.path.join(input_dir, sample_label)

        for f,filename in enumerate(os.listdir(subdirectory_path)):
            if ".txt" not in filename:
                continue
            # if "supported" in filename:
            #     continue
            if "nonref_node_length_vs" in filename:
                continue

            path = os.path.join(subdirectory_path, filename)
            series_label = filename.replace(".txt",'')

            values = [0]

            with open(path) as file:
                for l,line in enumerate(file):
                    x = float(line.strip())

                    values.append(x)

            print(series_label)
            print(values[:10], min(values), max(values))

            range = data_ranges[series_label]
            h,bins = numpy.histogram(values, 300, range)
            # h = h/100000
            h = h/h.sum()
            print(numpy.diff(bins))

            bins = numpy.cumsum(numpy.diff(bins))

            print(h)
            print(sum(h))
            print()

            data[sample_label][series_label] = (h,bins)

    # TODO: ratio of supported/unuspported vars by length

    # plot_6_major_distributions(data)
    plot_all(data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_dir",
        required=True,
        type=str,
        help="Path to directory containing the summary CSVs for a given evaluation"
    )

    args = parser.parse_args()

    print(args)

    main(args.input_dir)
