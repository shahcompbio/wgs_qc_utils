from classifycopynumber import parsers, transformations
import wgs_analysis.algorithms.cnv
import numpy as np



remixt_filename = snakemake.input['remixt_cn']
filename = snakemake.output[0]
sample = snakemake.wildcards[0]

cn, stats = parsers.read_remixt_parsed_csv(remixt_filename)
cn = cn.astype({"chromosome": "str"})

cn["sample"] = [sample] * len(cn)

aggregated_cn_data = wgs_analysis.algorithms.cnv.aggregate_adjacent(
    cn,
    value_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2'],
    stable_cols=['major_0', 'minor_0', 'major_1', 'minor_1', 'major_2', 'minor_2', 'sample'],
    length_normalized_cols=['major_raw', 'minor_raw'],
)
aggregated_cn_data['copy'] = aggregated_cn_data['major_raw'] + aggregated_cn_data['minor_raw']

aggregated_cn_data['ploidy'] = stats["ploidy"]
aggregated_cn_data['seg.mean'] = np.log2(aggregated_cn_data['copy'] 
    / aggregated_cn_data['ploidy']
)
aggregated_cn_data['num.mark'] = (aggregated_cn_data['length'] / 500000).astype(int)
aggregated_cn_data = aggregated_cn_data.rename(
    columns={'sample': 'ID', 'chromosome': 'chrom', 
        'start': 'loc.start', 'end': 'loc.end'
    }
)
aggregated_cn_data = aggregated_cn_data[['ID', 'chrom', 'loc.start', 
    'loc.end', 'num.mark', 'seg.mean'
]]
aggregated_cn_data['seg.mean'] = aggregated_cn_data['seg.mean'].fillna(np.exp(-8))
aggregated_cn_data.loc[aggregated_cn_data['seg.mean'] == np.NINF, 'seg.mean'] = np.exp(-8)
aggregated_cn_data = transformations._correct_seg_bin_ends(aggregated_cn_data)
aggregated_cn_data.to_csv(filename, index=None, sep='\t')

