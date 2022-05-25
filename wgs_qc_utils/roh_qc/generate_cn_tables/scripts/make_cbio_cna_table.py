import pandas as pd

cn_change_filename = snakemake.input[0]
cbio_table = snakemake.output[0]

chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']

gistic_data = pd.read_csv(cn_change_filename, converters={'chromosome': str})
gistic_data = gistic_data[gistic_data['chromosome'].isin(chromosomes)]
gistic_data = gistic_data.rename(columns={"gene_name": "Hugo_Symbol"})
gistic_data = gistic_data.astype({"gistic_value": "Int64"})

# Duplicate entries are possible, select the entry with the highest absolute gistic value
gistic_data['abs_gistic_value'] = gistic_data['gistic_value'].abs()
abs_max_gistic = gistic_data.groupby(['Hugo_Symbol', 'sample'])['abs_gistic_value'].max().reset_index()
gistic_data = gistic_data.merge(abs_max_gistic).drop_duplicates(['Hugo_Symbol', 'sample', 'gistic_value'])

gistic_matrix = gistic_data.set_index(['Hugo_Symbol', 'sample'])['gistic_value'].unstack()
gistic_matrix.reset_index(inplace=True)
gistic_matrix.to_csv(cbio_table, sep="\t", index=False, na_rep="NA")
