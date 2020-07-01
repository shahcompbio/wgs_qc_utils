import pandas as pd
import numpy as np


def read(coverage, binsize=100000):
    """
    read in coverage data
    :param coverage: coverage data
    :return: pandas dataframe
    """

    cov = pd.read_csv(coverage, na_values="nan", sep="\t", 
                      names=["chrom", "start", "end", "sum_cov"])
    

    cov = cov.astype({"chrom": str, "start": "int64", "end":"int64", "sum_cov":"int64"})

    cov["coverage"] = cov.sum_cov/binsize


    return cov

def prepare_at_chrom(coverage, chrom, bin=False, n_bins=200):
    """
    prepare coverage data at a chromosome to be plotted
    :param coverage: read-in coverage data
    :param chrom: chromosome
    :param n_bins: number of bins to bin data into
    :return: binned dataframe
    """
    coverage = coverage[coverage.chrom == str(chrom)]

    if bin:
        coverage = bin_coverages(coverage, n_bins)

    high = coverage.coverage.quantile(0.99)
    low = coverage.coverage.quantile(0.01)

    # preserve NaNs so as not to plot over centromere
    return coverage[(coverage.coverage.between(low, high))
                    | (coverage.coverage.isnull())]


def bin_coverages(coverage, n_bins):
    """
    bin coverage data
    :param coverage: input coverage data
    :param n_bins: number of bins to seperate data into
    :return: binned coverage data as dataframe
    """
    positions = coverage.start
    coverages = coverage.coverage

    bins = np.linspace(positions.min(), positions.max(), n_bins)
    digitized = np.digitize(positions, bins)

    start = [positions[digitized == i].min() for i in range(1, len(bins))]
    end = [positions[digitized == i].max() for i in range(1, len(bins))]
    coverage = [coverages[positions[digitized == i].index].mean() for i in range(1, len(bins))]

    return pd.DataFrame({"start": start, "end": end, "coverage": coverage})

