import pandas as pd
import snv_cn
import remixt_plotting
import variant_plotting
import pysam

def tri(bam, chrom, pos):
    bam.fetch(chrom, pos-1, pos+1)
    return bam.sequence

remixt_9 = remixt_plotting.read("/Users/abramsd/work/DATA/QC/remixt/results_files_009.h5", "009")
somatic_calls = variant_plotting.read_consensus("/Users/abramsd/work/DATA/QC/somatic/Sample009_somatic.csv.gz")
vaf_data = snv_cn.parse(somatic_calls, remixt_9)
vaf_data.to_csv("vaf009.csv", index=False)
chroms = vaf_data.chr.unique()
bam = pysam.AlignmentFile("")

for chrom in chroms:
    vaf = snv_cn.prepare_at_chrom(vaf_data, chrom)
    vaf["trin"] = pd.pos.apply(lambda pos: tri(bam, chrom, pos))

freqs = pd.value_counts(vaf.trin)
print(freqs)
print(vaf_data)
