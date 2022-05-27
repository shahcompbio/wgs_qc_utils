# ROH QC

Creates a single-page PDF including ROH status for all chromosomes.

Run the following command (test run version):

```
snakemake \
  --configfile config.yaml \
  --cores 12
```

You can directly edit current config for your `isabl_sample_id`:

```
snakemake \
  --configfile config.yaml \
  --cores 12 \
  --config isabl_sample_id=ISABL_SAMPLE_ID
```

TODO: Dockerize to enable something like:

```
CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

snakemake --jobs 500 --skip-script-cleanup \
  --configfile config.yaml \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /juno"
```

