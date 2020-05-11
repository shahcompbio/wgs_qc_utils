import pandas as pd
import vcf
import gzip
import numpy as np


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
    print(f)
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
       'DP_tum':"float64",  })
    
    data[["allele_1", "allele_2"]] = data.GT_norm.str.split("/", expand=True)

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


def parse(reader, sep):
    """
    parse vcf
    :param f: vcf file
    :param sep: seperator to use
    :return: parsed vcf file as list of lists
    """
    return [line.strip().split(sep) for line in reader
            if not line.startswith("#")]
