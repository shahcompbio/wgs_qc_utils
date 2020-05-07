import titan_plotting
import matplotlib.pyplot as plt
import variant_plotting
import ideogram_plotting
import roh_plotting
import gene_annotation_plotting

fig, axes = plt.subplots(3,1, gridspec_kw={'height_ratios': [1,1,0.25]})
titan_7 = titan_plotting.read("/Users/abramsd/work/DATA/QC/titan/Sample007_titan_markers.csv.gz")
# remixt_7 = remixt_plotting.read("/Users/abramsd/work/DATA/QC/remixt/results_files_007.h5", "007")
somatic_calls = variant_plotting.read_consensus("/Users/abramsd/work/DATA/QC/somatic/Sample_007_somatic.csv.gz")
ideo = ideogram_plotting.read("/Users/abramsd/Downloads/cytoBandIdeo.txt")
roh = roh_plotting.read("/Users/abramsd/work/DATA/QC/roh/007_roh_ST.txt")
# vaf_data = snv_cn.parse(somatic_calls, remixt_7)

chrom = "1"
anno = gene_annotation_plotting.get_gene_annotation_data(chrom)
titan = titan_plotting.prepare_at_chrom(titan_7, chrom)
# remixt = remixt_plotting.prepare_at_chrom(remixt_7, chrom)
# vaf_data = snv_cn.prepare_at_chrom(vaf_data, chrom)
ideo = ideogram_plotting.prepare_at_chrom(ideo, chrom)
roh = roh_plotting.prepare_at_chrom(roh, chrom)

print(titan, roh, ideo)
titan_plotting.plot(titan,  anno, axes[0], titan.Position.max()/1000000)
roh_plotting.plot(roh,  axes[1], titan.Position.max()/1000000)
ideogram_plotting.plot(ideo, axes[2])
# snv_cn.plot_scatter(vaf_data, axes[0][0], logistic_y=False)
# snv_cn.plot_hist(vaf_data, axes[0][1], logistic_y=False)
# remixt_plotting.plot(remixt,  axes[1][0], remixt.start.max()/1000000, logistic_y=True)
# snv_cn.plot_scatter(vaf_data, axes[1][0], logistic_y=True)
# snv_cn.plot_hist(vaf_data, axes[1][1], logistic_y=True)

plt.show()
# plt.savefig("logistic_y_remixt.png")