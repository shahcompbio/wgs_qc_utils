# Generate CN tables pipeline

Prior to kicking off the command, make sure to create a `logs` directory in the same directory as the snakemake file in order for the logs to be saved.
The `logs` directory can also be configured in the `config.yaml` file under `log_dir:` to any other location.

Run the following command:

```
CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}")

snakemake --jobs 500 --skip-script-cleanup \
  --cluster-config cluster.yaml \
  --cluster "${CLUSTER_CMD}" \
  --use-singularity \
  --singularity-args "--bind /juno"
```

