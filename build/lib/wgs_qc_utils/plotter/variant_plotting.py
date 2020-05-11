import pandas as pd
import numpy as np
from . import abs_checker
import seaborn
import matplotlib
import wgs_analysis.wgs_analysis.refgenome as refgenome
import wgs_analysis.wgs_analysis.plots.colors
import wgs_analysis.wgs_analysis.plots as plots
import wgs_analysis.wgs_analysis.annotation.position as position

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

def rainfall(chrom, pos, ax):
    abs_checker.check_input_is_valid([chrom, pos],
                                     [abs_checker.CheckerTypes.STRING,
                                      abs_checker.CheckerTypes.INT])
      
    snvs = pd.DataFrame({"chrom":chrom.tolist(), "coord":pos.tolist()})
    snvs = position.annotate_adjacent_distance(snvs)
    snvs = snvs.drop_duplicates(['chrom', 'coord'])

    snvs = snvs.loc[(snvs['chrom'].isin(refgenome.info.chromosomes))]

    snvs.set_index('chrom', inplace=True)
    snvs['chromosome_start'] = refgenome.info.chromosome_start
    snvs['chromosome_color'] = pd.Series(plots.colors.create_chromosome_color_map())
    snvs.reset_index(inplace=True)

    snvs['plot_coord'] = snvs['coord'] + snvs['chromosome_start']

    assert not snvs['chromosome_color'].isnull().any()

    snvs['adjacent_distance_log'] = snvs['adjacent_distance'].apply(np.log10)

    ax.scatter(
        snvs['plot_coord'],
        snvs['adjacent_distance_log'], 
        facecolors=list(snvs['chromosome_color']),
        edgecolors=list(snvs['chromosome_color']),
        alpha=0.5, s=5, lw=0)

    ax.set_xlabel('chromosome')
    ax.set_xlim(min(snvs['plot_coord']), max(snvs['plot_coord']))
    ax.set_xticks([0] + list(wgs_analysis.wgs_analysis.refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(wgs_analysis.wgs_analysis.refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(wgs_analysis.wgs_analysis.refgenome.info.chromosomes))

    ax.set_ylim(-0.0001, 1.1 * snvs['adjacent_distance_log'].fillna(0.).max())
    ax.set_ylabel('Distance between mutations (log10)')

    seaborn.despine(trim=True)
        
    return ax
    

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