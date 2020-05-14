import pandas as pd
from . import abs_checker
import wgs_analysis.plots.snv as snv_plots
import wgs_analysis.annotation.position as position


def plot_bar(location, n_events, axis, name, chrom_max):
    """
    bar plot of variants on axis
    :param location: location col from read in variant calls
    :param n_events: n_events col from read in variant calls
    :param axis: axis to plot on
    :param name: name of y axis label
    :param chrom_max: max of x axis
    :return:  updated axis
    """

    abs_checker.check_input_is_valid([location, n_events],
                                     [abs_checker.CheckerTypes.NUMERIC,
                                      abs_checker.CheckerTypes.INT])

    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    axis.bar(location / 1000000,  n_events,
             facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    return axis


def plot_fill(location, n_events, axis, name, chrom_max):
    '''
    filled line plot of variants on axis
    :param location: location col from read in variant calls
    :param n_events: n_events col from read in variant calls
    :param axis: axis to plot on
    :param name: name of y axis label
    :param chrom_max: max of x axis
    :return:  updated axis
    '''
    abs_checker.check_input_is_valid([location, n_events],
                                     [abs_checker.CheckerTypes.NUMERIC,
                                      abs_checker.CheckerTypes.INT])

    axis.set_xlim(0, chrom_max)

    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.set_xticklabels([])
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False


    axis.fill_between(location / 1000000, 0, n_events,
                      facecolor='black', alpha=0.5)

    axis.set_ylabel(name, fontsize=14, fontname="Arial")

    return axis


def plot_rainfall(chrom, pos, axis):
    """
    rainfall plot
    :param chrom:
    :param pos:
    :param axis:
    :return:
    """
    abs_checker.check_input_is_valid([chrom, pos],
                                     [abs_checker.CheckerTypes.STRING,
                                      abs_checker.CheckerTypes.INT])

    snvs = pd.DataFrame({"chrom":chrom.tolist(), "coord":pos.tolist()})
    snvs = position.annotate_adjacent_distance(snvs)

    snv_plots.snv_adjacent_distance_plot(axis, snvs)

    return axis


