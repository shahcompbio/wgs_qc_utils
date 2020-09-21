from . import input_checker


def plot(pos, state, axis, chrom_max):
    """
    plot roh data on axes
    :param roh: prepped roh data
    :param axis: axis to plot on
    :return: axis with plot
    """


    if state.empty:
        return axis
        
    input_checker.check_input_is_valid([pos, state],
                                            [input_checker.CheckerTypes.NUMERIC,
                                             input_checker.CheckerTypes.NUMERIC])

    axis.fill_between(pos / 1000000, state.min(), state,
                      facecolor='black', alpha=0.5)
    axis.set_ylabel("ROH", fontsize=14, fontname="Arial")
    axis.set_xlim(0, chrom_max)
    axis.set_ylim(0, 1)
    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    return axis
