#!/usr/bin/env python

import os
import sys
import re
import gzip
import vcf
import click
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator, FixedFormatter

import scgenome.refgenome as refgenome


def annotate_size_class(data):
    """ Add the categorical column size_class.
    Expects columns:
        - position_1
        - position_2
        - type
    """

    length = (data['position_1'] - data['position_2']).abs()
    data['size_class'] = pd.cut(
        length,
        [0, 1e4, 1e6, 1e8, 1e10],
        labels=['0-1K', '1K-1M', '1-100M', '>100M'])
    data['size_class'] = data['size_class'].astype(str)
    data.loc[data['type'] == 'translocation', 'size_class'] = 'Tr'
    data['size_class'] = pd.Categorical(
        data['size_class'],
        categories=reversed(['0-1K', '1K-1M', '1-100M', '>100M', 'Tr']))

    return data


def classify_rearrangement_type(data):
    """Add the categorical column rearragement_type.
    Expects columns:
        - position_1
        - position_2
        - type
    """
    assert 'position_1' in data
    assert 'position_2' in data
    assert 'type' in data
    
    data = data.copy()
    
    break1 = data['position_1']
    break2 = data['position_2']
    size = abs(break1 - break2)
    orientation_type = data['type']
    
    data['rearrangement_type'] = 'unbalanced'
    
    deletion_ix = ((size <= 1e6) & (orientation_type == 'deletion'))
    foldback_ix = (~deletion_ix) & ((size <= 1e3) & (orientation_type == 'inversion'))
    inversion_ix = (~foldback_ix) & ((size <= 1e6) & (orientation_type == 'inversion'))
    duplication_ix = (~inversion_ix) & ((size <= 1e6) & (orientation_type == 'duplication'))
    # unbalanced_ix = (~duplication_ix)
    
    data.loc[deletion_ix, 'rearrangement_type'] = 'deletion'
    data.loc[foldback_ix, 'rearrangement_type'] = 'foldback'
    data.loc[inversion_ix, 'rearrangement_type'] = 'inversion'
    data.loc[duplication_ix, 'rearrangement_type'] = 'duplication'
    # data.loc[unbalanced_ix, 'rearrangement_type'] = 'unbalanced'
    
    return data


def create_breakends(breakpoints, data_cols=(), id_col='grouped_breakpoint_id'):
    breakends = breakpoints[[id_col, 'chromosome_1', 'position_1', 'chromosome_2', 'position_2']].copy()
    breakends.set_index(id_col, inplace=True)
    breakends.columns = pd.MultiIndex.from_tuples([tuple(c.split('_')) for c in breakends.columns])
    breakends = breakends.stack()
    breakends.index.names = (id_col, 'prediction_side')
    breakends = breakends.reset_index()
    breakends['prediction_side'] = np.where(breakends['prediction_side'] == '1', 0, 1)
    breakends = breakends.merge(breakpoints[[id_col] + list(data_cols)], on=id_col)
    return breakends


def type_size_plot(ax, data, type_col='type', log_scale=True):
    data = data.copy()
    data['type'] = data['type'].astype('category')
    data['rearrangement_type'] = pd.Categorical(
        data['rearrangement_type'],
        categories=[
            'duplication',
            'inversion',
            'balanced',
            'foldback',
            'deletion',
            'unbalanced',
        ])
    data['type'] = pd.Categorical(
        data['type'],
        categories=[
            'duplication',
            'inversion',
            'deletion',
            'translocation',
        ])

    p = sns.countplot(ax=ax, data=data,
        x=type_col, hue='size_class')
    
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    if log_scale:
        p.set_yticklabels()
        p.set(yscale="log")


def chromosome_type_plot(ax, breakends, bin_size=20000000, rearrangement_types=None, 
                         full_genome=True):
    default_rearrangement_types = ['foldback', 'deletion', 'duplication', 'inversion', 'balanced', 'unbalanced', 'complex']
    default_rearrangement_type_colors = sns.color_palette('Dark2', len(default_rearrangement_types))
    
    if rearrangement_types is None:
        rearrangement_types = default_rearrangement_types
        rearrangement_type_colors = default_rearrangement_type_colors
    else:
        rearrangement_type_colors = sns.color_palette('Dark2', len(rearrangement_types))

    breakends = breakends.loc[(breakends['chromosome'].isin(refgenome.info.chromosomes))].copy()

    breakends.set_index('chromosome', inplace=True)
    breakends['chromosome_start'] = refgenome.info.chromosome_start
    breakends.reset_index(inplace=True)

    breakends['plot_position'] = breakends['position'] + breakends['chromosome_start']
    breakends['plot_bin'] = breakends['plot_position'] / bin_size
    breakends['plot_bin'] = breakends['plot_bin'].astype(int)

    if full_genome:
        plot_bin_starts = np.arange(0, refgenome.info.chromosome_end.max(), bin_size)
    else:
        min_x, max_x = breakends['plot_position'].min(), breakends['plot_position'].max()
        plot_bin_starts = np.arange(min_x, max_x, bin_size)
    plot_bins = (plot_bin_starts / bin_size).astype(int)

    breakends = breakends[['grouped_breakpoint_id', 'plot_bin', 'rearrangement_type']].drop_duplicates()

    count_table = breakends.groupby(['plot_bin', 'rearrangement_type']).size().unstack()
    count_table = count_table.reindex(plot_bins).fillna(0)

    accum_counts = None
    for idx, rearrangement_type in enumerate(rearrangement_types):
        if rearrangement_type not in count_table.columns:
            continue
        ax.bar(plot_bin_starts, count_table[rearrangement_type].values,
               bottom=accum_counts, width=bin_size,
               facecolor=rearrangement_type_colors[idx],
               edgecolor=rearrangement_type_colors[idx])
        if accum_counts is None:
            accum_counts = count_table[rearrangement_type].values
        else:
            accum_counts += count_table[rearrangement_type].values

    if full_genome:
        ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
    else:
        ax.set_xlim((min_x, max_x))

    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.set_ylabel('count')
    ax.xaxis.set_minor_locator(FixedLocator(refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(FixedFormatter(refgenome.info.chromosomes))
    ax.xaxis.grid(False, which="minor")
    chromosome_type_plot_legend(ax)


def chromosome_type_plot_legend(ax, rearrangement_types=None):
    default_rearrangement_types = ['foldback', 'deletion', 'duplication', 'inversion', 'balanced', 'unbalanced', 'complex']
    default_rearrangement_type_colors = sns.color_palette('Dark2', len(default_rearrangement_types))
    
    if rearrangement_types is None:
        rearrangement_types = default_rearrangement_types
        rearrangement_type_colors = default_rearrangement_type_colors
    else:
        rearrangement_type_colors = sns.color_palette('Dark2', len(rearrangement_types))

    ax.legend([plt.Circle((0, 0), color=c) for c in rearrangement_type_colors],
              list(rearrangement_types), loc="upper right")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))


def draw_per_size_plot(data, ax):
    # per-type plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    type_size_plot(ax, data, type_col='rearrangement_type', log_scale=False)


def draw_per_chrom_plot(data, ax):
    # per-chromosome plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    breakends = create_breakends(breakpoints=data, 
                                 id_col='grouped_breakpoint_id',
                                 data_cols=['rearrangement_type'])
    chromosome_type_plot(ax, breakends)
    

def draw_hist_and_ecdf(data, axes):
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    # Draw PDE
    splt = sns.histplot(data=data, ax=axes[0],
                        x='n_split', hue="rearrangement_type", 
                        element = 'poly', fill=False, 
                        binrange=(0, 50), bins=20)
    splt.set(xlim=(0, None))

    # Draw CDE
    splt = sns.ecdfplot(data=data, ax=axes[1],
                        x='n_split', hue="rearrangement_type")
    splt.set(xlim=(0, 50))


def classify_alt_svtype(alt):
    assert type(alt) == str
    alt_pos_regex = re.search('([0-9A-Z\.]+):([0-9]+)', alt)
    if alt == '<DEL>':
        svtype = 'deletion'
    elif alt == '<DUP>':
        svtype = 'duplication'
    elif alt == '<INV>':
        svtype = 'inversion'
    elif alt_pos_regex:
        svtype = 'translocation'
    else:
        svtype = None # error
    return svtype


def get_destruct_data(destruct_path):
    destruct = pd.read_table(destruct_path, sep='\t', converters={
        'chromosome_1': str,
        'chromosome_2': str,
    })
    chrom1s = destruct.chromosome_1.values
    chrom2s = destruct.chromosome_2.values
    pos1s = destruct.position_1.astype(int).values
    pos2s = destruct.position_2.astype(int).values
    ixs = zip(chrom1s, pos1s, chrom2s, pos2s)
    num_splits = destruct.num_unique_reads.values
    svtypes = destruct.type.values
    num_splits_svtypes = zip(num_splits, svtypes)
    destruct_data = dict(zip(ixs, (num_splits_svtypes)))
    return destruct_data


def get_svaba_data(svaba_path):
    svaba_data = {}
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 
              'FORMAT', 'V1', 'V2', 'V3', 'NORMAL', 'TUMOR']
    pos1s = []
    pos2s = []
    drs = []
    for line in gzip.open(svaba_path, 'rb'):
        line = line.decode("utf-8")
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            # header = line.strip().split('\t') # -> svaba output is full-o-shit
            continue
        field = line.strip().split('\t')
        row = dict(zip(header, field))
        chrom1 = row['#CHROM']
        pos1 = int(row['POS'])
        formats = row['FORMAT'].split(':')
        s1s = row['NORMAL'].split(':')
        s2s = row['TUMOR'].split(':')
        tumor = dict(zip(formats, s2s))
        dr = int(tumor['AD']) # Allele depth: Number of reads supporting the variant

        alt = str(row['ALT']) # assume 1 alt
        alt_pos_regex = re.search('([0-9A-Z\.]+):([0-9]+)', str(alt))
        if alt_pos_regex:
            chrom2, pos2 = alt_pos_regex.groups()
            pos2 = int(pos2)
        else:
            error_row = row
            print(f"ERROR: {row} does not have ending position")
            break

        svtype = classify_alt_svtype(alt)
        if svtype == None:
            print('[ERROR] alt =', alt)
            break

        svaba_data[(chrom1, pos1, chrom2, pos2)] = (dr, svtype)
    return svaba_data


def get_lumpy_data(lumpy_path, sample):
    lumpy_data = {}

    vcf_reader = vcf.Reader(filename=lumpy_path)
    for row in vcf_reader:
        sample_rows = [s for s in row.samples if s.sample == sample]
        assert len(sample_rows) == 1
        sample_row = sample_rows[0]
        chrom1 = str(row.CHROM)
        pos1 = int(row.POS)
        alt = str(row.ALT[0]) # assume 1 alt
        alt_pos_regex = re.search('([0-9A-Z\.]+):([0-9]+)', str(alt))
        if 'END' in row.INFO:
            chrom2 = str(row.CHROM)
            pos2 = int(row.INFO['END'])
        elif alt_pos_regex:
            chrom2, pos2 = alt_pos_regex.groups()
            pos2 = int(pos2)
        else:
            error_row = row
            print(f"ERROR: {row} does not have ending position")
            break

        svtype = classify_alt_svtype(alt)
        if svtype == None:
            print('[ERROR] alt =', alt)
            break

        su = sample_row['SU'] # Number of pieces of evidence supporting the variant
        lumpy_data[(chrom1, pos1, chrom2, pos2)] = (su, svtype)
    return lumpy_data


def get_gridss_data(gridss_path):
    gridss_data = {}

    vcf_reader = vcf.Reader(filename=gridss_path)
    for row in vcf_reader:
        if len(row.FILTER) > 0:
            continue
        sample_rows = [s for s in row.samples if s.sample == 'tumour']
        assert len(sample_rows) == 1
        sample_row = sample_rows[0]
        chrom1 = str(row.CHROM)
        pos1 = int(row.POS)
        alt = str(row.ALT[0]) # assume 1 alt
        alt_pos_regex = re.search('([0-9A-Z\.]+):([0-9]+)', str(alt))
        if alt_pos_regex:
            chrom2, pos2 = alt_pos_regex.groups()
            pos2 = int(pos2)
        else:
            continue # single breakends
        svtype = classify_alt_svtype(alt)
        if svtype == None:
            print('[ERROR] alt =', alt)
            break

        vf = sample_row['VF'] # Count of fragments supporting the variant breakpoint allele and not the reference allele
        gridss_data[(chrom1, pos1, chrom2, pos2)] = (vf, svtype)
    return gridss_data


def make_sv_table_from_data(sv_data):
    chromosome_1s = []
    position_1s = []
    chromosome_2s = []
    position_2s = []
    types = []
    n_splits = []
    grouped_breakpoint_ids = []
    id_ix = 0
    for key, value in sv_data.items(): # ('1', '85028338', '1', '85028329'): (1, 'inversion'), ...
        id_ix += 1
        grouped_breakpoint_ids.append(id_ix)
        chrom1, pos1, chrom2, pos2 = key
        n_split, svtype = value
        chromosome_1s.append(chrom1)
        position_1s.append(pos1)
        chromosome_2s.append(chrom2)
        position_2s.append(pos2)
        types.append(svtype)
        n_splits.append(int(n_split))

    sv = pd.DataFrame(zip(chromosome_1s, position_1s, 
                          chromosome_2s, position_2s,
                          types, n_splits, grouped_breakpoint_ids),
                      columns = ['chromosome_1', 'position_1', 
                                 'chromosome_2', 'position_2',
                                 'type', 'n_split', 'grouped_breakpoint_id'])
    return sv


def update_consensus_with_read_count(consensus, destruct_data, svaba_data, 
        lumpy_data, gridss_data):
    # Filter out variants with <N supporting reads
    consensus_copy = consensus.copy()
    n_splits = []
    destruct_svtype, svaba_svtype, lumpy_svtype, gridss_svtype = None, None, None, None
    for ix, row in consensus_copy.iterrows():
        n_split_aggregate = []
        key = (row.chromosome_1, row.position_1, 
               row.chromosome_2, row.position_2)
        if key in destruct_data:
            destruct_n_split, destruct_svtype = destruct_data[key]
            n_split_aggregate.append(destruct_n_split)
        if key in svaba_data:
            svaba_n_split, svaba_svtype = svaba_data[key]
            n_split_aggregate.append(svaba_n_split)
        if key in lumpy_data:
            lumpy_n_split, lumpy_svtype = lumpy_data[key]
            n_split_aggregate.append(lumpy_n_split)
        if key in gridss_data:
            gridss_n_split, gridss_svtype = gridss_data[key]
            n_split_aggregate.append(gridss_n_split)
        if destruct_svtype == svaba_svtype == lumpy_svtype == gridss_svtype == None:
            sys.exit('ERROR: all svtypes are None')
        n_splits.append(max(n_split_aggregate))
    consensus_copy['n_split'] = pd.Series(n_splits)
    return consensus_copy
# Get input


def make_plots_pdf(pdf, min_split_reads_cutoffs, svs, sample):
    with PdfPages(pdf) as pdf_pages:
        page_width = 12 # inches
        page_height = 9 # inches
        
        for tool in ('consensus', 'destruct', 'lumpy', 'svaba', 'gridss',):
            print(f'\n[LOG] processing tool: {tool}\n', file=sys.stderr)
            table = svs[tool]

            fig = plt.figure(figsize=(page_width, page_height))
            plt.axis('off')
            plt.title(f'\n\nBreakpoint QC ({tool}): {sample}\n\n', fontweight='bold')
            plt.subplots_adjust(hspace=0.5, top=0.7, bottom=0.3)
            gs = gridspec.GridSpec(nrows=1, ncols=2)
            ax00 = fig.add_subplot(gs[0, 0])
            ax01 = fig.add_subplot(gs[0, 1])
            ax00.set_ylabel("frequency")
            ax00.set_xlabel("split reads")
            ax01.set_ylabel("CDF")
            ax01.set_xlabel("split reads")
            
            for ix, min_split_reads_cutoff in enumerate(min_split_reads_cutoffs):
                print(f'[LOG] processing min_split_reads_cutoff: {min_split_reads_cutoff}', file=sys.stderr)
                flt_table = table[table.n_split >= min_split_reads_cutoff]
                data = annotate_size_class(flt_table)
                data = classify_rearrangement_type(data)

                if min_split_reads_cutoff == 0:
                    draw_hist_and_ecdf(data, axes=(ax00, ax01,))
                    pdf_pages.savefig(fig)
                    plt.close()

                fig = plt.figure(figsize=(page_width, page_height))
                plt.axis('off')
                plt.title(f"{tool}: breakpoints with supporting reads ≥{min_split_reads_cutoff}")
                plt.subplots_adjust(hspace=0.5, right=0.8)
                gs = gridspec.GridSpec(nrows=2, ncols=1)

                ax0 = fig.add_subplot(gs[0, 0])
                ax1 = fig.add_subplot(gs[1, 0]) 
                draw_per_size_plot(data, ax0) 
                draw_per_chrom_plot(data, ax1)
                pdf_pages.savefig(fig)
                plt.close()


@click.command()
@click.option("--in_dir", type=click.Path(exists=True), required=True,
        help="input files directory")
@click.option("--sample", required=True,
        help="Isabl sample ID (str)")
@click.option("--pdf", required=True,
        help="output plot pdf path")
@click.option("--genome_version", default='hg19', show_default=True,
        help="select genome version")
def breakpoint_qc(in_dir, sample, pdf, genome_version='hg19'):
    #genome_version = 'mm10'
    #sample = 'MPCdk12_CL_early' # TODO: argumentize
    #in_dir = f'/juno/work/shah/users/grewald/tickets/SHAH-3365/breakpoint/{sample}/results/' # TODO directory change

    refgenome.set_genome_version(genome_version)

    consensus_path = f'{in_dir}/four_way_consensus.csv.gz' # consensus calls
    destruct_path = f'{in_dir}/{sample}_breakpoint_table.csv'
    lumpy_path = f'{in_dir}/{sample}_lumpy.vcf'
    svaba_path = f'{in_dir}/{sample}.svaba.somatic.sv.vcf.gz'
    #gridss_path = f'{in_dir}/{sample}_gridss.vcf.gz'
    gridss_path = f'{in_dir}/gridss.pass.test.vcf'

    assert os.path.isfile(consensus_path)
    assert os.path.isfile(consensus_path)
    assert os.path.isfile(destruct_path)
    assert os.path.isfile(lumpy_path)
    assert os.path.isfile(svaba_path)

    consensus = pd.read_csv(consensus_path,
            converters = {'chromosome_1': str, 'chromosome_2': str,
                          'position_1': int, 'position_2': int})
    destruct_data = get_destruct_data(destruct_path)
    svaba_data = get_svaba_data(svaba_path)
    lumpy_data = get_lumpy_data(lumpy_path, sample)
    gridss_data = get_gridss_data(gridss_path)
    consensus = update_consensus_with_read_count(consensus, 
            destruct_data, svaba_data, lumpy_data, gridss_data)

    svs = {
        'consensus': consensus,
        'destruct': make_sv_table_from_data(destruct_data),
        'svaba': make_sv_table_from_data(svaba_data),
        'lumpy': make_sv_table_from_data(lumpy_data),
        'gridss': make_sv_table_from_data(gridss_data),
    }

    #pdf = './test.pdf' # TODO argumentize
    min_split_reads_cutoffs = (0, 5, 10) # TODO: argumentize
    make_plots_pdf(pdf, min_split_reads_cutoffs, svs, sample)


if __name__ == "__main__":
    breakpoint_qc()
