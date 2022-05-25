import pandas as pd

sample = snakemake.wildcards[0]

df = pd.read_csv(snakemake.input[0], dtype='str')
df['sample'] = sample
df.to_csv(str(snakemake.output[0]), index=False)
