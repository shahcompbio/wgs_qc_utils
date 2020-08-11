library(maftools);
library(R.utils);


#m = "/Users/meredith/DOUGAS_WORK/work/DATA/pseudobulkQC/spectrum_maf.maf.gz"
m = "/juno/work/shah/svatrt/cbioportal_tools/example/apollo_outputs/data_mutations_extended.maf"
#m = "/juno/work/shah/abramsd/CODE/out2.maf"
#m = "/juno/work/shah/abramsd/CODE/apollo_maf_filtgenes_patients.maf"
#m="/juno/work/shah/svatrt/cbioportal_tools/example/diljot_outputs/data_mutations_extended.maf"
#csv = read.csv(m, sep="\t")
#csv <- csv[,colSums(is.na(csv))<nrow(csv)]

#print(csv)

#laml = read.maf(maf = m);
#print(laml)
#print(colnames(laml@data))
##write.plotmafSummary(maf = laml, basename="output", rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#csv = read.csv(m, sep="\t")
#csv <- csv[,colSums(is.na(csv))<nrow(csv)]
#csv = csv[csv$IMPACT=="HIGH",]
#write.csv(csv, "TESTEST", sep="\t")
#m = "TESTEST"
laml = read.maf(maf = m);
#oncoplot(maf = laml, top=20)
#print(getSampleSummary(laml))

genelist = read.csv("/juno/work/shah/abramsd/CODE/apollooncoplotgenelist.tsv")
genes = unlist(list(genelist$gene_name))

pdf("apollo_oncoplot_filtgenes_patients_high_impact.pdf")

#onocplot(maf=laml)
#pdf("spectrum_default_100_oncoplot")
#oncoplot(maf = laml, genes=c("TP53", "BRCA1", "BRCA2", "PIK3CA",
                            # "RB1", "PTEN", "CDK12", "PALB2",
                            # "MYC", "KRAS", "CCNE1", "NF1"))
genes =c("TPF3", "BRCA1", "CDK12", "NF1", "BRCA2", "RHGAP26", "BTK", "KMT2C", "NF2", "PBRM1")

oncoplot(maf=laml, genes=genes)
dev.off()
