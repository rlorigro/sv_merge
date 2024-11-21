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
    test_1074= [
        "/home/ryan/data/test_hapestry/run/1074_skip_logs/output"
    ]

    test_local= [
        "/home/ryan/data/test_hapestry/run/test_logs"
    ]

    times = defaultdict(lambda: defaultdict(list))

    # aggregate_logs(parent_dirs=test_1074, times_by_name=times["test_1074"])
    aggregate_logs(parent_dirs=test_local, times_by_name=times["test_local"])

    fig = pyplot.figure()
    ax = pyplot.axes()

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

    colormap = matplotlib.colormaps["nipy_spectral"]

    shuffled_indexes = list(range(0,len(names)))
    # random.shuffle(shuffled_indexes)

    for s,solver in enumerate(times):
        bottom = 0

        for n,name in enumerate(names):
            y = sum(times[solver][name])

            y = float(y/60)
            x = s + 0.5

            l = str(n) + "_" + name

            color_index = (float(shuffled_indexes[n]) + 1) /float(len(names))

            print(n, len(names) , color_index)

            # color_index = float(n)/float(len(names)-1)
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

            if "total" not in name:
                x_t = (float(x) - w/2.0) + (color_index * w/1.0)
                y_t = (b + b + y) / 2
                ax.text(x_t, y_t, "-" + str(n) + "-",va='center', ha='center')

                x_ticks.append(x)
                x_labels.append(solver)

                bottom += y

    pyplot.xticks(x_ticks, x_labels, rotation=45, ha="right")

    ax.set_xlabel("Solver")
    ax.set_ylabel("Time (m)")

    pyplot.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    pyplot.tight_layout()

    pyplot.show()
    pyplot.close()



if __name__ == "__main__":
    main()
