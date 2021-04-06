import pandas as pd
import numpy as np
from wgs_qc_utils.reader.read_variant_calls import EmptyVariantReader

class EmptySnvCnReader():
    def __init__(self):
        self.pos = None
        self.frac_cn = None


def parse(snvs, remixt):
    if isinstance(snvs, EmptyVariantReader):
        return EmptySnvCnReader()
    snv_cn_table = annotate_copy_number(snvs, remixt,
                                        columns=['major', 'minor', 'total_raw_e',
                                                 'tumour_content', 'is_subclonal'])
    # print(snv_cn_table, snv_cn_table.columns)
    output =  pd.DataFrame([calculate_cellular_frequency(row) for i, row in snv_cn_table.iterrows()])
    output.rename(columns={"chromosome":"chrom"})
    return output


def prepare_at_chrom(transformed_snv, chrom):
    if isinstance(transformed_snv, EmptySnvCnReader):
        return EmptySnvCnReader()
    out = transformed_snv[transformed_snv.chrom == chrom]
    return out


def annotate_copy_number(pos, seg, columns=['major', 'minor']):
    results = []

    for chrom in seg['chrom'].unique():
        _pos = pos[pos['chrom'] == chrom]
        _seg = seg[seg['chrom'] == chrom]
        results.append(find_overlapping_segments(_pos, _seg, columns))
    return pd.concat(results)


def find_overlapping_segments(pos, seg, columns):
    seg = seg.sort_values(['start', 'end'])
    if seg.duplicated().any():
        raise ValueError('duplicate columns')
    start_idx = np.searchsorted(seg['start'].values, pos['pos'].values) - 1
    end_idx = np.searchsorted(seg['end'].values, pos['pos'].values)
    mask = (start_idx == end_idx)
    results = pos.copy()
    for col in columns:
        results[col] = np.nan
        results.loc[mask, col] = seg[col].iloc[end_idx[mask]].values
    return results


def calculate_cellular_frequency(row):
    if row['major'] == 0:
        row['ccf'] = 0.
        row['alt_cn'] = 0
        return row
    ccfs = list()
    cns = list()
    total_depth = (2. * (1. - row['tumour_content']) + row['total_raw_e'] * row['tumour_content'])
    for cn in range(1, int(row['major']) + 1):
        alt_depth = row['tumour_content'] * cn
        ccf = row['VAF_tumor'] * total_depth / alt_depth
        ccfs.append(ccf)
        cns.append(cn)
    ccfs = np.array(ccfs)
    cns = np.array(cns)
    idx = np.argmin(np.absolute(1. - ccfs))
    row['ccf'] = ccfs[idx]
    row['alt_cn'] = cns[idx]
    row['frac_cn'] = row['VAF_tumor'] * total_depth / row['tumour_content']
    return row
