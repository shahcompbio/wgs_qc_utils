import pandas as pd

data = pd.read_csv("/juno/work/shah/abramsd/oncokb-annotator/apollo_merged_maf_oncokb_anno.maf", sep="\t")
data["patient_id"]  = data.apply(lambda r: r.Tumor_Sample_Barcode.split("AD")[1][1:], axis=1)
data["Tumor_Sample_Barcode"] = data.patient_id
data.to_csv("/juno/work/shah/abramsd/oncokb-annotator/apollo_merged_maf_oncokb_patients.maf", sep="\t", index=False)


