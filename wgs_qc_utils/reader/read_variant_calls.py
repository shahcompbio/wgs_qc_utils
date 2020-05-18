import pandas as pd
import vcf
import gzip
import numpy as np


def prepare_at_chrom(variants, chrom, n_bins=200):
    '''
    prepare variants data to be plotted at a chrom
    '''

    variants = variants[variants["chr"] == str(chrom)]

    return bin_frequencies(variants.pos, n_bins, variants.pos.min(),
                           variants.pos.max())


def bin_frequencies(locations, n_bins, start, extent):
    '''
    bin variant data
    '''
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(locations, bins)

    binned_loc = [locations[digitized == i].mean() for i in range(1, len(bins))]
    n_events = [len(locations[digitized == i]) for i in range(1, len(bins))]

    return pd.DataFrame({"location": binned_loc,
                       "n_events": n_events})


def read_consensus_csv(f):
    f = pd.read_csv(f)
    f["VAF_normal"] = f.AC_NORMAL/f.DP_NORMAL
    f["VAF_tumor"] = f.AC_TUMOUR/f.DP_TUMOUR
    f = f.rename(columns = {"chrom":"chr"})
    f = f.astype({"chr":str})
    return f


def read_full_slow(f):
    reader = vcf.Reader(open(f))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out


def read_titan_vcf(f):
    if f.endswith(".gz"):
        lines = parse(gzip.open(f, "rt"), "\t")
    else:
        lines = parse(open(f), "\t")
    data = pd.DataFrame(lines, columns=["chr", "pos", "id", "ref", "alt",
                                                          "qual", "filter", "info",
                                                           "format", "tumour", "normal"])

    data = data.astype({"pos": np.int64, "chr": str})[["chr", "pos","format",
                                                           "tumour", "normal"]]
    cols = data.format.str.split(":").tolist()[0]

    cols_norm = [c + "_norm" for c in cols]
    data[cols_norm] = data.normal.str.split(":", expand=True)
    cols_tum= [c + "_tum" for c in cols]
    data[cols_tum] = data.tumour.str.split(":", expand=True)

    data = data.drop("format", axis=1)
    data = data.drop("normal", axis=1)

    data = data.astype({'RC_norm':"float64", 'AC_norm':"float64", 'NI_norm':"float64", 'ND_norm':"float64",
       'DP_norm':"float64",   'RC_tum':"float64", 'AC_tum':"float64", 'NI_tum':"float64", 'ND_tum':"float64",
       'DP_tum':"float64"})
    
    data[["allele_1", "allele_2"]] = data.GT_norm.str.split("/", expand=True)
    #
    data["normal_vaf"] = data.RC_norm/data.DP_norm
    data["tumour_vaf"] = data.RC_tum/data.DP_tum

    return data


def read(f):
    '''
    read in
    '''
    if f.endswith(".gz"):
        lines = parse(gzip.open(f, "rt"), "\t")
    else:
        lines = parse(open(f), "\t")

    data = pd.DataFrame(lines,
                        columns=["chr", "pos", "id", "ref", "alt", "qual",
                                 "filter", "info", "format", "normal"])

    data = data.astype({"pos": np.int64, "chr": str})

    cols = data.format.str.split(":").tolist()[0]
    data[cols] = data.normal.str.split(":", expand=True)

    data = data.drop("format", axis=1)
    data = data.drop("normal", axis=1)

    return data


def read_svs(breakpoints):
    breakpoints = pd.read_csv(breakpoints)[["chromosome_1", "chromosome_2", "position_1", "position_2", "prediction_id", "rearrangement_type"]]
    breakpoints = breakpoints.astype({"chromosome_1": str, "chromosome_2": str, "rearrangement_type":str})

    breakpoints = pd.DataFrame({"chr": breakpoints["chromosome_1"].append(breakpoints["chromosome_2"]),
                                "pos": breakpoints["position_1"].append(breakpoints["position_2"]),
                                "rearrangement_type": breakpoints["rearrangement_type"].append(breakpoints["rearrangement_type"]),
                                "prediction_id": breakpoints["prediction_id"].append(breakpoints["prediction_id"])})
    return breakpoints


def parse(reader, sep):
    """
    parse vcf
    :param f: vcf file
    :param sep: seperator to use
    :return: parsed vcf file as list of lists
    """
    return [line.strip().split(sep) for line in reader
            if not line.startswith("#")]

# variant_file = "/Users/abramsd/work/DATA/QC/titan/Sample007_museq.vcf"
# data = read_titan_vcf(variant_file)
# import matplotlib.pyplot as plt
# plt.hist(data.normal_vaf, bins = 100)
# print(data[data.normal_vaf > 100000][["normal_vaf", "DP_norm", "RC_norm"]])
# plt.show()