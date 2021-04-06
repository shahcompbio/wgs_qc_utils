import pandas as pd
from . import input_checker
import logging
import os
from wgs_qc_utils.utils.empty import empty_plot


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
    if not isinstance(location, pd.Series) and not isinstance(n_events, pd.Series):
        return empty_plot(axis, "Breakpoints")

    input_checker.check_input_is_valid([location, n_events],
                                     [input_checker.CheckerTypes.NUMERIC,
                                      input_checker.CheckerTypes.INT])

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


def plot_fill(location, n_events, axis, name, chrom_max, type):
    '''
    filled line plot of variants on axis
    :param location: location col from read in variant calls
    :param n_events: n_events col from read in variant calls
    :param axis: axis to plot on
    :param name: name of y axis label
    :param chrom_max: max of x axis
    :return:  updated axis
    '''
    if not isinstance(location, pd.Series) and not isinstance(n_events, pd.Series):
        return empty_plot(axis, type)
        
    input_checker.check_input_is_valid([location, n_events],
                                     [input_checker.CheckerTypes.NUMERIC,
                                      input_checker.CheckerTypes.INT])

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
