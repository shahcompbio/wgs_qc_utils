import pandas as pd
import numpy as np


def read(coverage, binsize=100000):
    """
    read in coverage data
    :param coverage: coverage data
    :return: pandas dataframe
    """
    with open(coverage, "r") as f:
        header  = f.readline()
    f.close()

    if "start" in header:
        cov = pd.read_csv(coverage, na_values="nan", sep="\t")
    else:
        cov = pd.read_csv(coverage, na_values="nan",names=['chrom', 'start', 'end', 'sum_cov'], sep="\t")

    cov = cov.astype({"chrom": str, "start": "int64", "end":"int64", "sum_cov":"int64"})

    cov["coverage"] = cov.sum_cov/binsize
    cov["chrom"] = cov.chrom.str.lower()
 
    return cov

def prepare_at_chrom(coverage, chrom):
    """
    prepare coverage data at a chromosome to be plotted
    :param coverage: read-in coverage data
    :param chrom: chromosome
    :param n_bins: number of bins to bin data into
    :return: binned dataframe
    """
    coverage = coverage[coverage.chrom == str(chrom)]


    high = coverage.coverage.quantile(0.99)
    low = coverage.coverage.quantile(0.01)

    # preserve NaNs so as not to plot over centromere
    return coverage[(coverage.coverage.between(low, high))
                    | (coverage.coverage.isnull())]
