from collections import defaultdict
import os

import matplotlib
from matplotlib import pyplot


def aggregate_logs(parent_dirs, times, common_windows, name):
    success_counter = [set(),set()]
    fail_reasons = defaultdict(int)

    windows = defaultdict(set)

    for directory in parent_dirs:
        for root, dirs, files in sorted(os.walk(directory)):
            for file in files:
                file_path = os.path.join(root, file)

                if "log.csv" in file_path:
                    window = root.strip('/').split("/")[-1]

                    with open(file_path, 'r') as file:
                        for l,line in enumerate(file):
                            if l == 0:
                                continue

                            task, h, m, s, ms, success, notes = line.strip().split(',')
                            windows[task].add(window)

                            if not bool(int(success)):
                                if task != "window_total":
                                    fail_reasons[task] += 1
                                    print(window, task, notes)

                                # continue

                            success_counter[int(success)].add(window)

                            s_total = float(h)*60*60 + float(m)*60 + float(s) + float(ms)/1000.0
                            times[name][task].append((window,s_total))

    print(name,  len(success_counter[1]), len(success_counter[0]), ','.join(["%s:%d"%(k,v) for k,v in fail_reasons.items()]))

    for task in windows:
        if len(common_windows[task]) > 0:
            common_windows[task] = common_windows[task].intersection(windows[task])

            # just for convenience of debugging or general edification
            # uncommon_windows = windows[task].difference(common_windows[task])
            # print(task, uncommon_windows)
        else:
            common_windows[task] = windows[task]

    return times, common_windows


def plot_barchart(times, task_names):
    x_ticks = list()
    x_labels = list()

    observed = set()

    colormap = matplotlib.colormaps["nipy_spectral"]

    fig = pyplot.figure()
    ax = pyplot.axes()

    for r,run_name in enumerate(times):
        bottom = 0

        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue

            # y = sum([t for w,t in times[run_name][task_name] if w in common_windows[task_name]])
            y = sum([t for w,t in times[run_name][task_name]])

            print('\t'.join(list(map(str,[run_name, task_name, y]))))

            x = r + 0.5

            l = str(t) + "_" + task_name

            color_index = (float(t) + 1) /float(len(task_names))

            color = colormap(color_index)

            b = bottom
            z = 1
            w = 0.8
            if "total" in task_name:
                b = 0
                z = 0
                color = "gray"
                w = 0.82

            if task_name not in observed:
                p = ax.bar(x, y, w, label=l, bottom=b, color=color, zorder=z)
                observed.add(task_name)
            else:
                p = ax.bar(x, y, w, bottom=b, color=color, zorder=z)

            if "total" not in task_name:
                x_t = (float(x) - w/2.0) + (color_index * w/1.0)
                y_t = (b + b + y) / 2
                ax.text(x_t, y_t, "-" + str(t) + "-",va='center', ha='center')

                x_ticks.append(x)
                x_labels.append(run_name)

                bottom += y

    pyplot.xticks(x_ticks, x_labels, rotation=45, ha="right")

    ax.set_xlabel("Solver")
    ax.set_ylabel("Time (m)")

    pyplot.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    pyplot.tight_layout()

    pyplot.show()
    pyplot.close()


def plot_scaling(times, task_names):
    x_ticks = list()
    x_labels = list()

    observed = set()

    colormap = matplotlib.colormaps["nipy_spectral"]

    fig = pyplot.figure()
    ax = pyplot.axes()

    y_per_task = defaultdict(list)
    x = list()

    ordering = lambda x: int(x[1:])

    # Assume runs are named n0001,n0002, etc
    for r,run_name in enumerate(sorted(times, key=ordering)):
        bottom = 0

        x.append(ordering(run_name))

        for t,task_name in enumerate(task_names):
            if task_name not in times[run_name]:
                continue

            # y = sum([t for w,t in times[run_name][task_name] if w in common_windows[task_name]])
            y = sum([t for w,t in times[run_name][task_name]])

            y_per_task[task_name].append(y)


    for i,(task_name,y) in enumerate(y_per_task.items()):
        color_index = (float(i) + 1) /float(len(task_names))

        color = colormap(color_index)

        z = 1
        w = 0.8
        if "total" in task_name:
            b = 0
            z = 0
            color = "gray"
            w = 0.82

        pyplot.plot(x,y,color=color, label=task_name, marker='o')

        x_ticks.append(x)

    pyplot.xticks(x)

    ax.set_xlabel("n samples")
    ax.set_ylabel("Time (m)")

    pyplot.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    pyplot.tight_layout()

    pyplot.show()
    pyplot.close()


def main():
    n0001 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0001/"
    ]
    n0002 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0002/"
    ]
    n0004 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0004/"
    ]
    n0008 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0008/"
    ]
    n0016 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0016/"
    ]
    n0032 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0032/"
    ]
    n0064 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0064/"
    ]
    n0128 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0128/"
    ]
    n0256 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0256/"
    ]
    n0512 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n0512/"
    ]
    n1024 = [
        "/home/ryan/data/test_hapestry/run/aou_scaling/n1024/"
    ]
    n0256_trgt = [
        "/home/ryan/data/test_hapestry/run/test_trgt_scaling"
    ]

    times = defaultdict(lambda: defaultdict(list))
    common_windows = defaultdict(set)

    # times, common_windows = aggregate_logs(parent_dirs=n0001, times=times, common_windows=common_windows, name="n0001")
    # times, common_windows = aggregate_logs(parent_dirs=n0002, times=times, common_windows=common_windows, name="n0002")
    # times, common_windows = aggregate_logs(parent_dirs=n0004, times=times, common_windows=common_windows, name="n0004")
    # times, common_windows = aggregate_logs(parent_dirs=n0008, times=times, common_windows=common_windows, name="n0008")
    # times, common_windows = aggregate_logs(parent_dirs=n0016, times=times, common_windows=common_windows, name="n0016")
    # times, common_windows = aggregate_logs(parent_dirs=n0032, times=times, common_windows=common_windows, name="n0032")
    # times, common_windows = aggregate_logs(parent_dirs=n0064, times=times, common_windows=common_windows, name="n0064")
    # times, common_windows = aggregate_logs(parent_dirs=n0128, times=times, common_windows=common_windows, name="n0128")
    # times, common_windows = aggregate_logs(parent_dirs=n0256, times=times, common_windows=common_windows, name="n0256")
    # times, common_windows = aggregate_logs(parent_dirs=n0512, times=times, common_windows=common_windows, name="n0512")
    # times, common_windows = aggregate_logs(parent_dirs=n1024, times=times, common_windows=common_windows, name="n1024")
    times, common_windows = aggregate_logs(parent_dirs=n0256_trgt, times=times, common_windows=common_windows, name="n0256_trgt")

    task_names = [
        # "variant_graph",
        "graphaligner",
        "align_reads_to_paths",
        # "feasibility_construct",
        # "feasibility_init",
        "feasibility",
        # "compress_transmap",
        # "optimize_d_prune_construct",
        # "optimize_d_prune_init",
        "optimize_d_prune",
        # "optimize_d_construct",
        # "optimize_d_init",
        # "optimize_d",
        # "optimize_d_prune_parse",
        # "optimize_d_plus_n_construct",
        # "optimize_d_plus_n_init",
        "optimize_d_plus_n",
        # "optimize_d_plus_n_parse",
        "window_total",
    ]

    plot_barchart(times, task_names)
    # plot_scaling(times, task_names)



if __name__ == "__main__":
    main()
