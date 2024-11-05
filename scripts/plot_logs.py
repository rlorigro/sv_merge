from collections import defaultdict
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

                if "log.csv" in file_path:
                    with open(file_path, 'r') as file:
                        for l,line in enumerate(file):
                            if l == 0:
                                continue

                            name, h, m, s, ms, success, notes = line.strip().split(',')

                            s_total = float(h)*60*60 + float(m)*60 + float(s) + float(ms)/1000.0
                            times_by_name[name].append(s_total)


def main():
    gurobi_dirs = [
        "/home/ryan/data/test_hapestry/run/dev_small_test_gurobi"
    ]

    gurobi_samplewise_dirs = [
        "/home/ryan/data/test_hapestry/run/dev_small_test_samplewise_gurobi"
    ]

    scip_dirs = [
        "/home/ryan/data/test_hapestry/run/dev_small_test_samplewise"
    ]

    scip_samplewise_dirs = [
        "/home/ryan/data/test_hapestry/run/dev_small_test"
    ]


    times = dict()
    times["gurobi"] = defaultdict(list)
    times["gurobi_samplewise"] = defaultdict(list)
    times["scip"] = defaultdict(list)
    times["scip_samplewise"] = defaultdict(list)

    aggregate_logs(parent_dirs=gurobi_dirs, times_by_name=times["gurobi"])
    aggregate_logs(parent_dirs=gurobi_samplewise_dirs, times_by_name=times["gurobi_samplewise"])
    aggregate_logs(parent_dirs=scip_dirs, times_by_name=times["scip"])
    aggregate_logs(parent_dirs=scip_samplewise_dirs, times_by_name=times["scip_samplewise"])

    fig = pyplot.figure()
    ax = pyplot.axes()


    colors = {
        "graphaligner": "#5F18A0",
        "feasibility": "#EDBF15",
        "n_min": "#698C00",
        "d_given_n": "#88B107",
        "d_given_n_min": "#88B107",
        "optimize_d": "#AFE213",
        "d_min": "#AFE213",
        "optimize_n_given_d": "#C2EB42",
        "optimize_n_given_d_min": "#C2EB42",
        "optimize_d_plus_n": "#D3F56B",
        "optimize_n_d_quadratic": "#D3F56B",
    }

    x_ticks = list()
    x_labels = list()

    for s,solver in enumerate(times):
        bottom = 0

        for n,name in enumerate(times[solver]):
            y = sum(times[solver][name])

            y = float(y/60)
            x = s + 0.5

            if s == 0:
                l = [name if name != "confident" else "graphaligner"]
                p = ax.bar(x, y, 0.8, label=l, bottom=bottom, color=colors[name])
            else:
                p = ax.bar(x, y, 0.8, bottom=bottom, color=colors[name])

            x_ticks.append(x)
            x_labels.append(solver)
            bottom += y

    pyplot.xticks(x_ticks, x_labels)

    ax.set_xlabel("Solver")
    ax.set_ylabel("Time (m)")

    pyplot.legend()

    pyplot.show()
    pyplot.close()



if __name__ == "__main__":
    main()
