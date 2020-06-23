import matplotlib
matplotlib.use('Agg')
from matplotlib.lines import Line2D
import numpy as np
from . import abs_checker


def plot(start, major_raw, minor_raw, axis, chrom_max, logistic_y=False):
    """
    plot remixt data on axis
    :param start: start column from parsed remixt h5
    :param major_raw: major_raw col from parsed remixt h5
    :param minor_raw: minor_raw col from parsed remit h5
    :param anno_genes: gene annotations to add
    :param axis: axis to plot on
    :return: axis with plot added
    """

    abs_checker.check_input_is_valid([start, major_raw, minor_raw],
                                            [abs_checker.CheckerTypes.NUMERIC,
                                             abs_checker.CheckerTypes.FLOAT,
                                             abs_checker.CheckerTypes.FLOAT])

    if logistic_y:
        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)
        maj = squash_f(major_raw) * 3
        min = squash_f(minor_raw) * 3
        yticks = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        axis.set_yticks(yticks_squashed)
        axis.set_yticklabels(ytick_labels)
        axis.set_ylim((-0.01, 1.01))
        axis.spines['left'].set_bounds(0, 1)

    else:
        axis.set_ylim(0, 8)
        axis.set_yticks(range(0, 9))
        maj = major_raw
        min = minor_raw

    axis.set_xlim(0, chrom_max)

    axis.grid(True, linestyle=':')
    axis.spines['left'].set_position(('outward', 5))
    axis.spines['bottom'].set_position(('outward', 5))
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticks(np.arange(0, start.max() / 1000000, 25))

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.scatter(start / 1000000,
                 maj,  c="red", s=20, marker="o", alpha=0.4)
    axis.scatter(start / 1000000,
                 min,  c="blue", s=20, marker="o", alpha=0.4)

    axis.set_ylabel("Remixt", fontsize=14, fontname="Arial")

    return axis


def add_remixt_legend(axis):
    '''
    add a remixt legend to an axis
    :param axis: matplotlib axis to add legend to
    :return: modified axis
    '''
    labels = ["major axis", "minor axis", "snv copy number"]

    lines = [Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='r', markersize=10, alpha=0.4),
             Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='b', markersize=10, alpha=0.4),
             Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='black', markersize=5, alpha=1)]

    axis.legend(lines, labels, ncol=1, loc="center left", title="Remixt",
                frameon=False, borderpad=0, borderaxespad=0)

    return axis
