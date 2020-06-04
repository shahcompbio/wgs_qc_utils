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

    return read



def prepare_at_chrom(copy_number, chrom):
    """
    get copy number data at a chromosome
    :param copy_number: read in copy number data
    :param chrom: chromosome to extract (str)
    :return: copy number data at chrom
    """
    '''
    prep copy number rdata for plotting at a chrom
    '''

    if isinstance(chrom, int):
        chrom = str(chrom)

    return copy_number[copy_number["Chr"] == chrom][["Position", "LogRatio", "color"]]


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
    bin the data at each chromosome, get rid of nans + unused cols
    :param copy_number: copy number input (titan_markers.csv.gz)
    :param outfile: output path for prepped data
    :return: prepped cn data ready to go for circos plot
    """

    cn = read(copy_number)

    chroms = cn.Chr.unique()
    output = dict(zip(chroms, [[]] * len(chroms)))

    for chrom in chroms:
        cn_at_chrom = cn[cn.Chr == chrom]
        out = bin_copy_number(cn_at_chrom, n_bins=200)

        out["Chr"] = [chrom] * len(out.index)
        output[chrom] = out

    output = pd.concat(output)
    output = output[~output.Position.isna()]
    print(outfile)
    cc
    output.to_csv(outfile, index=False, header=True, sep="\t")