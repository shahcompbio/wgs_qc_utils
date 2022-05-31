# Breakpoint QC

Creates a multi-page PDF including breakpoint QC output.
Specifically, creates QC plot outputs per-tool (consensus, destruct, svaba, lumpy, gridss), and per-supporting-read-cutoff.

TODO list:
- allow arguments (currently a sample Mondrian mouse run is hard-coded)
- test out on a Mondrian result
- appify

Run the following command:

```
CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

snakemake --jobs 500 --skip-script-cleanup \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /juno"
```

