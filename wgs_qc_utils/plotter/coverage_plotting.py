from . import abs_checker


def plot(start, coverage, ylim_min, ylim_max, axis, name, chrom_max):
    """
    plot coverage data on an axis
    :param prepped_coverage: prepped coverage data
    :param ylim_min: min for y axis
    :param ylim_max: max for y axis
    :param axis: axis to plot on
    :param name: name for axis (i.e. normal or tumour coverage)
    :return: axis with plot
    # """

    abs_checker.check_input_is_valid([start, coverage],
                                     [abs_checker.CheckerTypes.NUMERIC,
                                      abs_checker.CheckerTypes.FLOAT])

    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.fill_between(start / 1000000, ylim_min,
                      coverage, facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    axis.set_ylim(ylim_min, ylim_max)

    return axis
