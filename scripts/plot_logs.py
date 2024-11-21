from collections import defaultdict
import random
import tarfile
import sys
import os

import matplotlib
from matplotlib import pyplot


def aggregate_logs(parent_dirs, times_by_name):
    for directory in parent_dirs:
        for root, dirs, files in os.walk(directory):
            for file in files:
                file_path = os.path.join(root, file)
                print(file_path)

                if "log.csv" in file_path and not "eval" in file_path:
                    with open(file_path, 'r') as file:
                        for l,line in enumerate(file):
                            if l == 0:
                                continue

                            name, h, m, s, ms, success, notes = line.strip().split(',')

                            s_total = float(h)*60*60 + float(m)*60 + float(s) + float(ms)/1000.0
                            times_by_name[name].append(s_total)


def main():
    # with_csvs = [
    #     "/home/ryan/data/test_hapestry/run/1074_skip_logs/with_csv"
    # ]
    #
    # skip_csvs = [
    #     "/home/ryan/data/test_hapestry/run/1074_skip_logs/skip_csv"
    # ]

    test_local = [
        "/home/ryan/data/test_hapestry/run/test_logs"
    ]

    times = defaultdict(lambda: defaultdict(list))

    aggregate_logs(parent_dirs=test_local, times_by_name=times["test_local"])
    # aggregate_logs(parent_dirs=skip_csvs, times_by_name=times["skip_csvs"])

    fig = pyplot.figure()
    ax = pyplot.axes()


    # colors = {
    #     "graphaligner": "#5F18A0",
    #     "variant_graph": "#AF0F87",
    #     "rescale_weights": "#19529C",
    #     "write_to_hap_vcf": "#0D9B77",
    #     "write_to_vcf": "#4FB69C",
    #     "align_reads_to_paths": "#3222A5",
    #     "feasibility": "#EDBF15",
    #     "n_min": "#698C00",
    #     "d_given_n": "#88B107",
    #     "d_given_n_min": "#88B107",
    #     "optimize_d": "#AFE213",
    #     "d_min": "#AFE213",
    #     "optimize_n_given_d": "#C2EB42",
    #     "optimize_n_given_d_min": "#D3F56B",
    #     "optimize_d_plus_n": "#C2EB42",
    #     "optimize_n_d_quadratic": "#D3F56B",
    #     "optimize_d_initial": "red",
    #     "window_total": "gray"
    # }

    names = [
        # "variant_graph",
        "graphaligner",
        "align_reads_to_paths",
        "feasibility_construct",
        "feasibility_init",
        "feasibility",
        "optimize_d_prune_construct",
        "optimize_d_prune_init",
        "optimize_d_prune",
        "optimize_d_prune_parse",
        "optimize_d_plus_n_construct",
        "optimize_d_plus_n_init",
        "optimize_d_plus_n",
        "optimize_d_plus_n_parse",
        "window_total",
    ]

    x_ticks = list()
    x_labels = list()

    observed = set()

    colormap = matplotlib.colormaps["hsv"]

    shuffled_indexes = list(range(0,len(names)))
    # random.shuffle(shuffled_indexes)

    for s,solver in enumerate(times):
        bottom = 0

        for n,name in enumerate(names):
            y = sum(times[solver][name])

            y = float(y/60)
            x = s + 0.5

            l = [name if name != "confident" else "graphaligner"]

            color_index = float(shuffled_indexes[n])/float(len(names))
            color = colormap(color_index)

            b = bottom
            z = 1
            w = 0.8
            if "total" in name:
                b = 0
                z = 0
                color = "gray"
                w = 0.82

            if name not in observed:
                p = ax.bar(x, y, w, label=l, bottom=b, color=color, zorder=z)
                observed.add(name)
            else:
                p = ax.bar(x, y, w, bottom=b, color=color, zorder=z)

            x_ticks.append(x)
            x_labels.append(solver)

            if "total" not in name:
                bottom += y

    pyplot.xticks(x_ticks, x_labels, rotation=45, ha="right")
    pyplot.tight_layout()

    ax.set_xlabel("Solver")
    ax.set_ylabel("Time (m)")

    pyplot.legend()

    pyplot.show()
    pyplot.close()



if __name__ == "__main__":
    main()
