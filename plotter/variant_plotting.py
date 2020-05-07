import pandas as pd
import numpy as np
import abs_checker

def plot_bar(location, n_events, axis, name, chrom_max):
    '''
    plot variants on axis
    '''

    abs_checker.check_input_is_valid([location, n_events],
                                     [abs_checker.CheckerTypes.INT,
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
    plot variants on axis
    '''

    abs_checker.check_input_is_valid([location, n_events],
                                     [abs_checker.CheckerTypes.INT,
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


def bin_frequencies(locations, n_bins, start, extent):
    '''
    bin variant data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(locations, bins)

    binned_loc = [locations[digitized == i].mean() for i in range(1, len(bins))]
    n_events = [len(locations[digitized == i]) for i in range(1, len(bins))]

    return pd.DataFrame({"location": binned_loc,
                       "n_events": n_events})


def prepare_at_chrom(variants, chrom, n_bins=200):
    '''
    prepare variants data to be plotted at a chrom
    '''

    variants = variants[variants["chr"] == str(chrom)]

    return bin_frequencies(variants.pos, n_bins, variants.pos.min(),
                           variants.pos.max())