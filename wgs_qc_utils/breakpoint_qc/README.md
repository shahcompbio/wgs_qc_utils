# Breakpoint QC

Creates a multi-page PDF including breakpoint QC output.
Specifically, creates QC plot outputs per-tool (consensus, destruct, svaba, lumpy, gridss), and per-supporting-read-cutoff.

TODO list:
- test out on a Mondrian result
- appify
- add example/test

Run the following command:

```
bsub -n 1 -J CLUSTER_JOB_NAME -o CLUSTER_LOG.out -e CLUSTER_LOG.err \
module load singularity/3.6.2 && \
singularity run docker://shahlab/breakpoint_qc:0.0.1 \
/path/to/breakpoint_qc.py \
  --in_dir /path/to/input/directory \
  --sample ISABL_SAMPLE_ID \
  --pdf /path/to/output.pdf
```

