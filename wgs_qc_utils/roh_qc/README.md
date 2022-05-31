# ROH QC

Creates a single-page PDF including ROH status for all chromosomes.

Run the following commands. You can directly edit current config for your `isabl_sample_id`, as in the example.

```
CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

module load singularity/3.6.2 && \
snakemake --jobs 500 --skip-script-cleanup \
  --configfile config.yaml \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /juno"
  --config isabl_sample_id=ISABL_SAMPLE_ID
```

