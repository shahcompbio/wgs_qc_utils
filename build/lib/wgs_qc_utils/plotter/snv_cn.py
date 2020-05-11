import pandas as pd
import numpy as np


def annotate_copy_number(pos, seg, columns=['major', 'minor']):
    results = []

    for chrom in seg['chromosome'].unique():
        _pos = pos[pos['chr'] == chrom]
        _seg = seg[seg['chromosome'] == chrom]
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


def parse(snvs, remixt):
    snv_cn_table = annotate_copy_number(snvs, remixt,
                                        columns=['major', 'minor', 'total_raw_e',
                                                 'tumour_content', 'is_subclonal'])
    return pd.DataFrame([calculate_cellular_frequency(row)
                         for i, row in snv_cn_table.iterrows()])


def prepare_at_chrom(transformed_snv, chrom):
    out = transformed_snv[transformed_snv.chr == chrom]
    out["frac_cn"] = out.frac_cn *3
    return out

def plot_scatter(data, axis, logistic_y=False):
    if logistic_y:
        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)
        frac_cn = squash_f(data.frac_cn)
    else:
        frac_cn = data.frac_cn
    axis.scatter(data.pos/1000000, frac_cn, color="black", s=2.5, marker="o", alpha=1)
    return axis


def plot_hist(data, axis, logistic_y=False):
    if logistic_y:
        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)
        frac_cn = squash_f(data.frac_cn)
        yticks = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        axis.set_yticks(yticks_squashed)
        axis.set_yticklabels(ytick_labels)
        axis.set_ylim((-0.01, 1.01))
        axis.spines['left'].set_bounds(0, 1)
    else:
        axis.set_ylim(0, 8)
        frac_cn = data.frac_cn
    axis.hist(frac_cn, color="black", orientation="horizontal", bins=15)
    axis.set_ylabel("SNV density")
    return axis
