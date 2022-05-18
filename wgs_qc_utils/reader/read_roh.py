import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d


class EmptyRohReader():
    def __init__(self):
        self.start = None
        self.state = None


def read(file):
    """
    read in roh data to pandas df
    :param file: roh file
    :return: pandas dataframe
    """
    if file.endswith(".gz"):
        data = pd.read_csv(
            file, 
            usecols=["type","sample","chromosome","start","state","quality"],
            converters={
                'chromosome': str,
            }
        )
    else:
        data = _parse_old_roh_format(file)
    if data.empty:
        return EmptyRohReader()
    data = data.rename(columns={"chromosome": "chrom", "quality":"qual"})
    data = data[~data.state.isna()]
    data = data.astype({"chrom":str})
    data["chrom"] = data.chrom.str.lower()
    
    return data
    

def _parse_old_roh_format(roh_calls): 
    lines = [l.strip("\n").split("\t") for l in open(roh_calls) if "ST" in l and "#" not in l]
    data = pd.DataFrame(lines,
        columns = ["type", "sample", "chromosome", "start", "state", "quality"]
    )
    data = data.astype({"type": "str", "sample": "str", "chromosome": "str",
        "start": "int", "state": "str", "quality": "float"}
    )
    return data


def prepare_at_chrom(roh, chrom):
    """
    prep roh data for plotting
    :param roh: read-in roh data
    :param chrom: chromosome
    :return: pandas dataframe
    """
    if isinstance(roh, EmptyRohReader):
        return roh
        
    roh = roh[(roh.qual > 30) & (roh.chrom == chrom)]
    smoothed_state = gaussian_filter1d(list(map(float, roh.state)), sigma=1000)
    smoothed_state = [0 if v < 0.49 else 1 for v in smoothed_state]
    roh["state"] = smoothed_state
    return roh
