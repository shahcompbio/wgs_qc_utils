import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d


def read(file):
    """
    read in roh data to pandas df
    :param file: roh file
    :return: pandas dataframe
    """
    data = pd.read_csv(file, usecols=["type","sample","chromosome","start","state","quality"])
    data = data.rename(columns={"chromosome": "chrom", "quality":"qual"})
    data = data.astype({"chrom":str})
    data["chrom"] = data.chrom.str.lower()
    
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
