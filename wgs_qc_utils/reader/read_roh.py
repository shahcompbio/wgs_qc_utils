import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d


def read(file):
    """
    read in roh data to pandas df
    :param file: roh file
    :return: pandas dataframe
    """

    return rename(pd.read_csv(file, sep="\t"))


def rename(data):
    """
    simplify column names
    :param data: read-in roh data dataframe
    :return: pandas dataframe
    """

    data = data.rename(columns={" # ST": "ST",
                                "[2]Sample": "sample",
                                "[3]Chromosome": "chrom",
                                "[4]Position": "pos",
                                "[5]State (0:HW, 1:AZ)": "state",
                                "[6]Quality (fwd-bwd phred score)": "qual"})
    data = data.astype({"chrom":str})
    return data


def prepare_at_chrom(roh, chrom):
    """
    prep roh data for plotting
    :param roh: read-in roh data
    :param chrom: chromosome
    :return: pandas dataframe
    """

    roh = roh[(roh.qual > 30) & (roh.chrom == chrom)]
    smoothed_state = gaussian_filter1d(list(map(float, roh.state)), sigma=1000)
    smoothed_state = [0 if v < 0.49 else 1 for v in smoothed_state]
    roh["state"] = smoothed_state
    return roh


#might want it in future
# def bin_roh(data, bins):
#     """
#     bin roh data
#     :param data: read-in roh data
#     :param bins: number of bins to bin into
#     :return: binned pandas dataframe
#     """
#
#     pos = data.pos
#     state = data.state
#
#     bins = np.linspace(data.pos.min(), data.pos.max(), bins)
#     assigned = np.digitize(pos, bins)
#     bin_range = range(1, len(bins))
#
#     bin_start = [pos[assigned == i].min() for i in bin_range]
#
#     bin_end = [pos[assigned == i].max() for i in bin_range]
#     bin_mean = [pos[assigned == i].mean() for i in bin_range]
#
#     # calculate the number of assigned hom states over total number and round it
#     bin_state = [round(len(state[assigned == i][state[assigned == i] == 1])
#                  / (len(state[assigned == i]) + 1)) for i in bin_range]
#
#     n_homo_locs = [len(state[assigned == i][state[assigned == i] == 1])
#                    for i in bin_range]
#
#     n_het_locs = [len(state[assigned == i][state[assigned == i] == 0])
#                    for i in bin_range]
#
#     return pd.DataFrame({"bin_pos_mean": bin_mean, "bin_pos_start": bin_start,
#                          "bin_pos_end": bin_end, "bin_overall_state": bin_state,
#                          "n_homo_locs": n_homo_locs, "n_het_locs": n_het_locs})
