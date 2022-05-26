import click
import logging
import numpy as np
import matplotlib
import pandas as pd
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
from matplotlib.collections import BrokenBarHCollection

from wgs_qc_utils.plotter import roh_plotting

from wgs_qc_utils.reader import read_roh

from wgs_qc_utils.reader.ideogram import read_ideogram



def ideogram_plot(ideogram, axis):
    
    xranges = ideogram[['start', 'width']].values
    colors = ideogram["color"].values
    collection = BrokenBarHCollection(xranges, (0, 1), facecolors=colors, alpha=0.1)
    axis.add_collection(collection)
    axis.set_xlim(0, ideogram.start.max())
    axis.get_yaxis().set_visible(False)
    axis.set_xticks(np.arange(0, ideogram.start.max(), 25))

    return axis


def roh_plot(pos, state, axis, chrom_max):
    """
    plot roh data on axes
    :param roh: prepped roh data (pandas.DataFrame)
    :param axis: axis to plot on (pyplot.axis)
    :return: axis with plot (pyplot.axis)
    """
#     if not isinstance(pos, pd.Series) and not isinstance(state, pd.Series):
#         return empty_plot(axis, "roh")
        
#     input_checker.check_input_is_valid([pos, state],
#                                             [input_checker.CheckerTypes.NUMERIC,
#                                              input_checker.CheckerTypes.NUMERIC])

    margin = 0.1
    axis.fill_between(pos / 1000000, 
                      state.min()+state*margin, 
                      state-state*margin,
                      facecolor='red', alpha=1)
    axis.set_xlim(0, chrom_max)
    axis.set_ylim(0, 1)
    for tic in axis.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    return axis


@click.command()
@click.option("--roh", type=click.Path(exists=True), required=True,
        help="processed bcftools roh output (csv.gz)")
@click.option("--sample", required=True,
        help="Isabl sample ID (str)")
@click.option("--pdf", required=True,
        help="output plot pdf path")
@click.option("--chromosomes", 
        default = ','.join([str(_) for _ in range(1,23)] + ['x', 'y']),
        show_default=True,
        help="string of chromosomes[,] (str)")
def plot_roh_on_ideogram(roh, chromosomes, sample, pdf):
    chromosomes = chromosomes.split(',')
    roh = read_roh.read(roh)
    logging.debug(f'roh.head() = \n', roh.head())
    ideogram = read_ideogram.read()
    logging.debug(f'ideogram.head() = \n', ideogram.head())
    width_ratios = [ideogram[ideogram.chrom==chrom.lower()].start.max() for chrom in chromosomes]
    width_ratios = width_ratios / max(width_ratios)

    with matplotlib.backends.backend_pdf.PdfPages(pdf) as pdf:

        fig = plt.figure(figsize=(15, 9))
        fig.suptitle(f'\n\nROH for sample ID: {sample}', fontsize=16)
        nrows = 8
        ncols = 3
        gs = fig.add_gridspec(nrows, ncols)
        axes = [plt.subplot(cell) for cell in gs]

        chromosomes = [chrom.lower() for chrom in chromosomes]

        for ix, chrom in enumerate(chromosomes):
            logging.info(f'chromosome {chrom} being plotted')
            prepped_roh = read_roh.prepare_at_chrom(roh, chrom)
            prepped_ideogram = read_ideogram.prepare_at_chrom(ideogram, chrom)
            chrom_max = prepped_ideogram.start.max()
            
            width_ratio = width_ratios[ix]
            pos = axes[ix].get_position()
            fixed_pos = [pos.x0, pos.y0, pos.width*width_ratio, pos.height*0.5]
            axes[ix].set_position(fixed_pos)
            
            axes[ix] = ideogram_plot(prepped_ideogram, axes[ix])
            axes[ix] = roh_plot(prepped_roh.start, prepped_roh.state, axes[ix], chrom_max)
            axes[ix].set_title(f'chr{chrom.upper()}')
            axes[ix].set_xlabel('   Mbp', fontsize=8, ha='left')
            axes[ix].xaxis.set_label_coords(1, -0.3)
            axes[ix].set_rasterized(True)
                
        # plt.tight_layout() # fuck tight layout - will screw up plt.Axes.set_position
        pdf.savefig(fig)


if __name__ == "__main__":
    plot_roh_on_ideogram()
