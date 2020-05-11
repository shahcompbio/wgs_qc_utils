import pandas as pd
import numpy as np
from matplotlib.lines import Line2D



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
