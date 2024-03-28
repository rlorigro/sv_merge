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


def _invert(x, limits):
    """inverts a value x on a scale from
    limits[0] to limits[1]"""
    return limits[1] - (x - limits[0])


def _scale_data(data, ranges):
    """scales data[1:] to ranges[0],
    inverts if the scale is reversed"""
    # for d, (y1, y2) in zip(data[1:], ranges[1:]):
    for d, (y1, y2) in zip(data, ranges):
        assert (y1 <= d <= y2) or (y2 <= d <= y1)

    x1, x2 = ranges[0]
    d = data[0]

    if x1 > x2:
        d = _invert(d, (x1, x2))
        x1, x2 = x2, x1

    sdata = [d]

    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        if y1 > y2:
            d = _invert(d, (y1, y2))
            y1, y2 = y2, y1

        sdata.append((d-y1) / (y2-y1) * (x2 - x1) + x1)

    return sdata


def set_rgrids(self, radii, labels=None, angle=None, fmt=None,
               **kwargs):
    """
    Set the radial locations and labels of the *r* grids.
    The labels will appear at radial distances *radii* at the
    given *angle* in degrees.
    *labels*, if not None, is a ``len(radii)`` list of strings of the
    labels to use at each radius.
    If *labels* is None, the built-in formatter will be used.
    Return value is a list of tuples (*line*, *label*), where
    *line* is :class:`~matplotlib.lines.Line2D` instances and the
    *label* is :class:`~matplotlib.text.Text` instances.
    kwargs are optional text properties for the labels:
    %(Text)s
    ACCEPTS: sequence of floats
    """
    # Make sure we take into account unitized data
    radii = self.convert_xunits(radii)
    radii = np.asarray(radii)
    rmin = radii.min()
    # if rmin <= 0:
    #     raise ValueError('radial grids must be strictly positive')

    self.set_yticks(radii)
    if labels is not None:
        self.set_yticklabels(labels)
    elif fmt is not None:
        self.yaxis.set_major_formatter(FormatStrFormatter(fmt))
    if angle is None:
        angle = self.get_rlabel_position()

    self.set_rlabel_position(angle)

    for t in self.yaxis.get_ticklabels():
        t.update(kwargs)

    return self.yaxis.get_gridlines(), self.yaxis.get_ticklabels()


# taken from:
# https://datascience.stackexchange.com/questions/6084/how-do-i-create-a-complex-radar-chart
class ComplexRadar():
    def __init__(self, fig, variables, ranges, n_ordinate_levels=6):
        angles = np.arange(0, 360, 360./len(variables))

        axes = [fig.add_axes([0.2,0.2,0.6,0.6],polar=True, label="axes{}".format(i)) for i in range(len(variables))]

        l, text = axes[0].set_thetagrids(angles,labels=variables)

        for i, ax in enumerate(axes):
            grid = np.linspace(*ranges[i], num=n_ordinate_levels)
            gridlabel = ["{}".format(round(x,4)) for x in grid]

            # if ranges[i][0] > ranges[i][1]:
            #     grid = grid[::-1] # hack to invert grid

            gridlabel[0] = "" # clean up origin

            set_rgrids(ax, grid, labels=gridlabel, angle=angles[i])
            ax.set_ylim(*ranges[i])

        for label, angle in zip(axes[0].get_xticklabels(), angles):
            x,y = label.get_position()
            lab = axes[0].text(x,y-0.1, label.get_text(), transform=label.get_transform(),ha=label.get_ha(), va=label.get_va())

            a = angle-90
            if angle == 270:
                a = 0
            lab.set_rotation(a)

        axes[0].set_xticklabels([])

        for ax in axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)

        # variables for plotting
        self.angle = np.deg2rad(np.r_[angles, angles[0]])
        self.ranges = ranges
        self.ax = axes[0]

    def plot(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.plot(self.angle, np.r_[sdata, sdata[0]], *args, **kw)

    def fill(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


def plot_6_major_distributions(histograms: dict, title, colors):
    axes_locations = {
        "alignment_identity_avg": (0,0),
        "haplotype_coverage_avg": (0,1),
        # "cluster_alignment_identity_avg": (0,0),
        # "cluster_coverage_avg": (0,1),
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

            h = numpy.append(numpy.array([0]), h)
            bins = numpy.append(numpy.array([0]), bins)

            cdf = numpy.cumsum(h)

            if is_reverse[series_label]:
                cdf = 1 - cdf

            c = "red"
            if sample_label not in colors:
                sys.stderr.write("WARNING: could not find color for %s, using default of red" % sample_label)
            else:
                c = colors[sample_label]

            axes[i][j].plot(bins,cdf, label=sample_label, alpha=0.7, linewidth=2, color=c)
            axes[i][j].set_xlabel(series_label)

            axes[i][j].set_ylim([0,1.05])

            if is_reverse[series_label] and s == 0:
                axes[i][j].invert_xaxis()

    axes[0][0].legend()

    figure.suptitle(title)

    pyplot.show()
    pyplot.close()


def plot_all(histograms: dict, colors):
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

            h = numpy.append(numpy.array([0]), h)
            bins = numpy.append(numpy.array([0]), bins)

            cdf = numpy.cumsum(h)

            if is_reverse[series_label]:
                cdf = 1 - cdf

            all_series.append(cdf)

            p = axes.plot(bins,cdf, label=sample_label, alpha=0.7, linewidth=0.8, color=colors[sample_label])

            # c = p[-1].get_color()
            # axes.plot(bins,h, alpha=0.7, linewidth=0.8, color=c)

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


def load_dir(input_dir: str, data: defaultdict[lambda: defaultdict[list]]):
    """
    :param input_dir:
    :param data: will be populated as data[sample_label][series_label] = (h,bins)
    :return: data updated with more values
    """

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

            if series_label not in data[sample_label]:
                data[sample_label][series_label] = [0]

            with open(path) as file:
                for l,line in enumerate(file):
                    x = 0

                    if line.strip() == "-nan":
                        continue
                    else:
                        x = float(line.strip())

                    if x < 0:
                        continue

                    data[sample_label][series_label].append(x)

    return data


def evaluate(input_dir):
    # c6,c3,c1,c8,c9,c0,c4,
    colors = {
        'bcftools': "C0",
        'truvari': "C9",
        'svmerger': "C2",
        'jasmine': "C1",
        'svimmer': "C3",
        'sniffles_joint': "C6",
        'dipcall': "C4",
        'reference': "gray",
    }

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
    data = load_dir(input_dir=input_dir, data=data)

    histograms = defaultdict(dict)

    for sample_label in data.keys():
        for series_label,values in data[sample_label].items():
            range = data_ranges[series_label]
            h,bins = numpy.histogram(values, 300, range)
            h = h/h.sum()

            bins = numpy.cumsum(numpy.diff(bins))

            histograms[sample_label][series_label] = (h,bins)

    # TODO: ratio of supported/unsupported vars by length
    plot_6_major_distributions(data, title=os.path.basename(input_dir), colors=colors)


def plot_radar(data, data_ranges, title, colors):
    radar_series_labels = [
        "alignment_identity_avg",
        "haplotype_coverage_avg",
        # "n_alignments",
        "nonref_nodes_fully_covered",
        "nonref_edges_covered",
        # "nonref_nodes",
    ]

    for key,range in data_ranges.items():
        print(key)
        print(range)

        range[0] = max(0,range[0]*0.9999)
        # range[1] = min(1,range[1]*1.001)

        print(range)
        print()

        data_ranges[key] = range

    radar_ranges = [data_ranges[x] for x in radar_series_labels]

    fig = pyplot.figure()
    radar = ComplexRadar(fig, radar_series_labels, radar_ranges)

    for s,sample in enumerate(data):
        y = [data[sample][k] for k in radar_series_labels]

        print(sample)
        print(radar_series_labels)
        print(y)
        radar.plot(y, label=sample, linewidth=3.2, alpha=0.7, color=colors[sample])
        # radar.fill(y, alpha=0.05, color=colors[sample])

    # axes.legend()
    radar.ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.))
    fig.suptitle(title)

    fig.set_size_inches(8,8)
    pyplot.savefig("radar.png", dpi=300)
    pyplot.show()
    pyplot.close()


def evaluate_subdirs(parent_dir: str, use_radar: bool):
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

    # c6,c3,c1,c8,c9,c0,c4,
    colors = {
        'bcftools': "C0",
        'truvari': "C9",
        'svmerger': "C2",
        'jasmine': "C1",
        'svimmer': "C3",
        'sniffles_joint': "C6",
        'dipcall': "C4",
        'reference': "gray",
    }

    # data = defaultdict(lambda: defaultdict(list))
    # histograms = defaultdict(dict)

    data_per_stratification = defaultdict(lambda: defaultdict(dict))
    histograms_per_stratification = defaultdict(lambda: defaultdict(dict))
    averages_per_stratification = defaultdict(lambda: defaultdict(dict))

    '''
    analysis_small_overlap/               <--- parent dir
        ├── chr10_analysis_small          <--- subdir
        │    ├── all_windows              <--- stratification label
        │    │    ├── bcftools
        │    │    │    ├── alignment_identity_avg.txt
        │    │    │    ├── anomalous_windows.bed
        │    │    │    ├── cluster_alignment_identity_avg.txt
        │    │    │    ├── cluster_coverage_avg.txt
        │    │    │    ├── haplotype_coverage_avg.txt
        │    │    │    ├── n_alignments.txt
    '''

    for subdir in os.listdir(parent_dir):
        subdir_path = os.path.join(parent_dir, subdir)

        if not os.path.isdir(subdir_path):
            continue

        for stratification_label in os.listdir(subdir_path):
            stratification_dir = os.path.join(subdir_path, stratification_label)

            print("Loading data: ", subdir, stratification_label)
            data_per_stratification[stratification_label] = load_dir(
                input_dir=stratification_dir,
                data=data_per_stratification[stratification_label]
            )

    if not use_radar:
        for stratification_label,data in data_per_stratification.items():
            for sample_label in data.keys():
                for series_label,values in data[sample_label].items():
                    range = data_ranges[series_label]
                    h,bins = numpy.histogram(values, 300, range)
                    h = h/h.sum()

                    bins = numpy.cumsum(numpy.diff(bins))

                    histograms_per_stratification[stratification_label][sample_label][series_label] = (h,bins)

            print(histograms_per_stratification[stratification_label].keys())
            print("Plotting: ", stratification_label)

            plot_6_major_distributions(
                histograms_per_stratification[stratification_label],
                title=stratification_label,
                colors=colors
            )

    else:
        for stratification_label,data in data_per_stratification.items():
            empirical_data_ranges = defaultdict(lambda: [1e9,-1e9])

            for sample_label in data.keys():
                for series_label,values in data[sample_label].items():
                    avg = numpy.mean(values)

                    averages_per_stratification[stratification_label][sample_label][series_label] = avg

                    if avg < empirical_data_ranges[series_label][0]:
                        empirical_data_ranges[series_label][0] = avg

                    if avg > empirical_data_ranges[series_label][1]:
                        empirical_data_ranges[series_label][1] = avg

            print(histograms_per_stratification[stratification_label].keys())
            print("Plotting radar: ", stratification_label)

            plot_radar(
                data=averages_per_stratification[stratification_label],
                data_ranges=empirical_data_ranges,
                title=stratification_label,
                colors=colors
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i","--input_dir",
        required=False,
        type=str,
        help="Path to directory containing the summary CSVs for a given evaluation"
    )

    parser.add_argument(
        "-p","--parent_dir",
        required=False,
        type=str,
        help="Path to directory containing many parallel input directories, all of which are to be evaluated"
    )

    parser.add_argument(
        "-r","--radar",
        action=argparse.BooleanOptionalAction,
        help="If invoked, produce radar plots on this data"
    )

    args = parser.parse_args()

    has_input_dir = len(args.input_dir) != 0 if args.input_dir is not None else False
    parent_dir = len(args.parent_dir) != 0 if args.parent_dir is not None else False

    if has_input_dir == parent_dir == 0:
        exit("ERROR: must provide one of input_dir or parent_dir but not both")

    if has_input_dir:
        evaluate(args.input_dir)
    else:
        evaluate_subdirs(args.parent_dir, args.radar)

