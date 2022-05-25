import dask.dataframe as dd


df = dd.read_csv(snakemake.input, dtype='str')
df.to_csv(str(snakemake.output[0]), index=False, single_file=True)
