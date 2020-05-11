#Douglas Abrams

import pandas as pd
import numpy as np

def sort_chromosome_names(chromosomes):
    def get_chromosome_key(chromosome):
        try:
            return (0, int(chromosome))
        except ValueError:
            return (1, chromosome)
    return [chromosome for chromosome in sorted(chromosomes, key=get_chromosome_key)]


def parse_remixt_data(cnv, minlength=1000, major_col='major_raw', minor_col='minor_raw',
                    chromosome=None, start=None, end=None):

    if chromosome is None and (start is not None or end is not None):
        raise ValueError('start and end require chromosome arg')

    # Ensure we dont modify the calling function's table
    cnv = cnv[['chromosome', 'start', 'end', major_col, minor_col]].copy()

    if 'length' not in cnv:
        cnv['length'] = cnv['end'] - cnv['start']

    # Restrict segments to those plotted
    if chromosome is not None:
        cnv = cnv[cnv['chromosome'] == chromosome]
    if start is not None:
        cnv = cnv[cnv['end'] > start]
    if end is not None:
        cnv = cnv[cnv['start'] < end]

    # Create chromosome info table
    chromosomes = sort_chromosome_names(cnv['chromosome'].unique())
    chromosome_length = cnv.groupby('chromosome')['end'].max()
    chromosome_info = pd.DataFrame({'length': chromosome_length}, index=chromosomes)

    # Calculate start and end in plot
    chromosome_info['end'] = np.cumsum(chromosome_info['length'])
    chromosome_info['start'] = chromosome_info['end'] - chromosome_info['length']
    chromosome_info['mid'] = (chromosome_info['start'] + chromosome_info['end']) / 2.

    if minlength is not None:
        cnv = cnv[cnv['length'] >= minlength]

    cnv.set_index('chromosome', inplace=True)
    cnv['chromosome_start'] = chromosome_info['start']
    cnv.reset_index(inplace=True)

    cnv = cnv[["chromosome","start", "major_raw", "minor_raw"]]
    return cnv

def prep_remixt(remixt_file, sample_label):
    cn_data = {}
    with pd.HDFStore(remixt_file) as store:
        stats = store['stats'].sort_values('elbo').iloc[-1]
        stats['sample'] = sample_label
        init_id = stats['init_id']
        cn_data[sample_label] = store[f'/solutions/solution_{init_id}/cn']

    for sample, cn in cn_data.items():
        plot_data = cn.query('length > 100000').query('minor_readcount > 100')

        plot_data = plot_data.astype({"chromosome": str})

        prepped = parse_remixt_data(plot_data)
        return prepped
        # prepped.to_csv(prepped_name, sep="\t", index=False)



#use this to turn h5 output from remixt into circos input
def test():
    #example input (twins data)

    filenames = {
        "007":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_007.h5",
        "009":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_009.h5",
        "026":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_026.h5",
        "031":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_031.h5",
        "036":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_036.h5",
        "044":"/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/results_files_044.h5"
    }
    for sample, filename in filenames.items():
        filename = filenames[sample]

        output = "/Users/abramsd/work/CODE/QC/WWWW/circos_remixt/prepped_remixt_"+sample+".tsv"
        prep_remixt(filename, sample, output)

# test()