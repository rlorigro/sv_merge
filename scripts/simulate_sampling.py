from scipy.stats import poisson
from matplotlib import pyplot
import matplotlib
from matplotlib.ticker import FuncFormatter
import numpy


def coverage_plot():
    sample_size = 10_000.0

    min_coverage = 0.0
    max_coverage = 1.0
    c_step_size = 0.01

    # Allele frequency
    min_frequency = 0.0
    max_frequency = 0.1
    f_step_size = 0.001

    # Coverage threshold which we arbitraily define as "bad"
    threshold = 100.0

    af_steps = int(round((max_frequency - min_frequency) / f_step_size))
    coverage_steps = int(round((max_coverage - min_coverage) / c_step_size))

    matrix = numpy.zeros([af_steps, coverage_steps])

    f_tick_labels = list()
    c_tick_labels = list()

    for j in range(0,coverage_steps):
        c = min_coverage + j*f_step_size
        c_tick_labels.append(str('%.3f' % c))

    for i in range(0,af_steps):
        f = min_frequency + i*f_step_size
        f_tick_labels.append(str('%.3f' % f))

        for j in range(0,coverage_steps):
            c = min_coverage + j*c_step_size

            matrix[i,j] = sample_size * f * c

    color_map = matplotlib.colormaps.get_cmap('viridis')
    color_map.set_over("white")
    color_map.set_under("white")

    pyplot.imshow(matrix, vmax=threshold, cmap=color_map, origin='lower')
    axes = pyplot.gca()

    axes.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: '%.3f' % (min_frequency + float(x)*f_step_size)))
    axes.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: '%.3f' % (min_coverage + float(x)*c_step_size)))

    cbar = pyplot.colorbar()
    cbar.ax.set_ylabel('Effective allele coverage', rotation=270)
    # pyplot.tight_layout()

    axes.set_ylabel("Allele frequency")
    axes.set_xlabel("Coverage per sample")

    pyplot.show()


def poisson_plot():
    sample_size = 10_000.0
    genome_length = 3_000_000_000.0
    read_length = 20_000.0
    sv_length = 10_000.0
    o = 1_000.0
    q = read_length - o

    min_coverage = 0.0
    max_coverage = 1.0
    c_step_size = 0.01

    # Allele frequency
    min_frequency = 0.0
    max_frequency = 0.1
    f_step_size = 0.001

    # Coverage threshold which we arbitraily define
    max_threshold = 1.0 - 1e-12
    min_threshold = 0.0 + 1e-12

    af_steps = int(round((max_frequency - min_frequency) / f_step_size))
    coverage_steps = int(round((max_coverage - min_coverage) / c_step_size))

    matrix = numpy.zeros([af_steps, coverage_steps])

    f_tick_labels = list()
    c_tick_labels = list()

    for j in range(0,coverage_steps):
        c = min_coverage + j*f_step_size
        c_tick_labels.append(str('%.3f' % c))

    for i in range(0,af_steps):
        f = min_frequency + i*c_step_size
        f_tick_labels.append(str('%.3f' % f))

        for j in range(0,coverage_steps):
            c = min_coverage + j*f_step_size

            n_reads = c * genome_length * sample_size / read_length
            mu = n_reads/(genome_length*sample_size) * (sv_length+q-o)*f*sample_size

            p = poisson.cdf(20, mu)

            matrix[i,j] = 1 - p

    color_map = matplotlib.colormaps.get_cmap('viridis')
    color_map.set_over("white")
    color_map.set_under("white")

    pyplot.imshow(matrix,
                  vmin=min_threshold,
                  vmax=max_threshold,
                  cmap=color_map,
                  origin='lower')

    axes = pyplot.gca()

    axes.get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: '%.3f' % (min_frequency + float(x)*f_step_size)))
    axes.get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: '%.3f' % (min_coverage + float(x)*c_step_size)))

    cbar = pyplot.colorbar()
    cbar.ax.set_ylabel('Effective allele coverage', rotation=270)
    # pyplot.tight_layout()

    axes.set_ylabel("Allele frequency")
    axes.set_xlabel("Coverage per sample")

    pyplot.show()


if __name__ == "__main__":
    # coverage_plot()
    poisson_plot()
