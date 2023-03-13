import pandas as pd
import numpy as np


def read(copy_number):
    """
    read in copy number data
    :param copy_number: copy_number input data (titan_markers.csv.gz)
    :returns: read-in copynumber pandas dataframe with plot colors
    """

    read = pd.read_csv(copy_number, sep="\t")

    
    read = read.astype({"Chr": str})

    n_extra = read.TITANstate.max() - 7

    colors = ["#00FF00", "#006400", "#0000FF", "#880000",
             "#BB0000", "#CC0000", "#DD0000", "#EE0000"] + ["#FF0000"] * n_extra

    read["color"] = read.TITANstate.apply(lambda state: colors[state])
    read["Chr"] = (read.Chr
            .str.replace('chr', '')
            .str.lower())
    read.rename(columns={"Chr":"Chrom"}, inplace=True)

    return read



def prepare_at_chrom(copy_number, Chrom):
    """
    get copy number data at a Chromosome
    :param copy_number: read in copy number data
    :param Chrom: Chromosome to extract (str)
    :return: copy number data at Chrom
    """
    '''
    prep copy number rdata for plotting at a Chrom
    '''

    return copy_number[copy_number["Chrom"] == Chrom][["Position", "LogRatio", "color"]]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ for circos plotting ~~~~~~~~~~~~~~~~~~~ #


def bin_copy_number(copy_number, n_bins=200):
    """
    bin copy number data
    :param positions:
    :param copy_number:
    :param state:
    :param n_bins:
    :param start:
    :param extent:
    :return:
    """
    positions = copy_number.Position
    state = copy_number.TITANstate
    log_ratio = copy_number.LogRatio

    bins = np.linspace(copy_number.Position.min(), copy_number.Position.max(), n_bins)
    digitized = np.digitize(positions, bins)

    position = [positions[digitized == i].mean() for i in range(1, len(bins))]
    lr = [log_ratio[positions[digitized == i].index].mean() for i in range(1, len(bins))]
    state = [state[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame({"Position": position, "LogRatio": lr,
                         "state": state})


def make_for_circos(copy_number, outfile):
    """
    prepare cn data for circos plot
    bin the data at each Chromosome, get rid of nans + unused cols
    :param copy_number: copy number input (titan_markers.csv.gz)
    :param outfile: output path for prepped data
    :return: prepped cn data ready to go for circos plot
    """

    cn = read(copy_number)

    Chroms = cn.Chrom.unique()
    output = dict(zip(Chroms, [[]] * len(Chroms)))

    for Chrom in Chroms:
        cn_at_Chrom = cn[cn.Chrom == Chrom]
        out = bin_copy_number(cn_at_Chrom, n_bins=200)

        out["Chrom"] = [Chrom] * len(out.index)
        output[Chrom] = out

    output = pd.concat(output)
    output = output[~output.Position.isna()]
    output = output.rename(columns = {"Chrom": "Chr"})
    output.to_csv(outfile, index=False, header=True, sep="\t")

