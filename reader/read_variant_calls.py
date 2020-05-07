import pandas as pd
import vcf
import gzip
import numpy as np


def read_consensus_csv(f):
    f = pd.read_csv(f)
    f = f[["chrom", "pos", "AC_NORMAL", "DP_NORMAL", "AC_TUMOUR", "DP_TUMOUR"]]
    f["VAF_normal"] = f.AC_NORMAL/f.DP_NORMAL
    f["VAF_tumor"] = f.AC_TUMOUR/f.DP_TUMOUR
    f = f.drop(["AC_NORMAL", "DP_NORMAL", "AC_TUMOUR", "DP_TUMOUR"], axis=1)
    f = f.rename(columns = {"chrom":"chr"})
    f = f.astype({"chr":str})
    return f


def read_full_slow(f):
    reader = vcf.Reader(open(f))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out


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
