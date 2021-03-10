import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d


def read(file):
    """
    read in roh data to pandas df
    :param file: roh file
    :return: pandas dataframe
    """
    data = pd.read_csv(file, names=["ST","sample","chrom","pos","end","state","length","num_markers","qual"])
    data = data.astype({"chrom":str})
    data["chrom"] = data.chrom.str.lower()
    
    return data


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
