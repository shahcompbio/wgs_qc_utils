import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.lines import Line2D
from . import gene_annotation_plotting
from . import abs_checker


def plot(position, log_ratio, color, axis, chrom_max, anno_genes=[]):
    """
    plot prepped copy number data on axis
    :param prepped_copy_number: prepped copy number data (read->prepare_at_chrom->
    :param anno_genes:
    :param axis:
    :return:
    """
    abs_checker.check_input_is_valid([position, log_ratio, color],
                                            [abs_checker.CheckerTypes.NUMERIC,
                                             abs_checker.CheckerTypes.FLOAT,
                                             abs_checker.CheckerTypes.STRING])
    axis.set_ylim(-4, 6)
    axis.set_xlim(0, chrom_max)
    axis.grid(True, linestyle=':')
    axis.spines['left'].set_position(('outward', 5))
    axis.spines['bottom'].set_position(('outward', 5))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_yticks(np.arange(-4, 6, 1))
    axis.set_xticks(np.arange(0, position.max() / 1000000, 25))

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.scatter(position / 1000000,
                 log_ratio, s=0.1, color=color)

    axis.set_ylabel("Titan", fontsize=14, fontname="Arial")

    if any(anno_genes):
        axis = gene_annotation_plotting.plot_anno_genes(anno_genes, *axis.get_ylim(), axis)

    return axis


def add_titan_legend(axis):
    labels = ["1", "2", "3", "4", "5",
              "6", "7", "8", "9+"]

    colors = ["#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", "#FF0000"]

    lines = [Line2D([0], [0], marker='o', color='w', label='major axis',
                    markerfacecolor=c, markersize=10) for c, _ in zip(colors, labels)]

    axis.legend(lines, labels, ncol=3, columnspacing=0.02,
                loc="center left", title="Titan", frameon=False,
                borderpad=0, borderaxespad=0)

    return axis
