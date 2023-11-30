from matplotlib import pyplot


def main():
    # path = "/home/ryan/data/test_hapestry/hts_debug/hts_net_c_2000_iter_querys_loop.txt"
    # path = "/home/ryan/data/test_hapestry/hts_debug/hts_net_c_2000_iter_regions.txt"
    path = "/home/ryan/data/test_hapestry/hts_debug/hts_net_chr20.txt"

    x = list()
    sent = list()
    received = list()

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            if len(line.strip()) == 0:
                continue

            line = line.strip().split('\t')

            if len(line) != 3:
                continue

            print(float(l)*0.001, line)
            x.append(float(l)*0.001)
            sent.append(float(line[1]))
            received.append(float(line[2]))

    fig = pyplot.figure()
    axes = pyplot.axes()
    axes.plot(x,sent, label="sent")
    axes.plot(x,received, label="received")

    axes.legend()
    axes.set_xlabel("Time (s)")
    axes.set_ylabel("Rate (kB/s)")

    axes.set_ylim([0,140000])

    fig.set_size_inches(12,8)
    pyplot.savefig("hts_net_test_whole_chr20.png")

    # axes.axvline(4.685)
    # axes.axvline(4.685+3.661)

    pyplot.show()


if __name__ == "__main__":
    main()
