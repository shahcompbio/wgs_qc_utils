# ROH QC

Creates a single-page PDF including ROH status for all chromosomes.

Run the following command:

```
CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

snakemake --jobs 500 --skip-script-cleanup \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /juno"
```

