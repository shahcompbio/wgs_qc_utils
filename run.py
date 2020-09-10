import sv_tool_consensus as SVc
import pandas as pd
lumpy =  "/juno/work/shah/tantalus/SC-3803/results/breakpoint_calling/sample_SA1256PP/lumpy_breakpoints.csv.gz"
svaba = "/juno/work/shah/users/grewald/SVBENCH/SVABA/tempdir.svaba.somatic.sv.vcf.gz"
destruct = "/juno/work/shah/tantalus/SC-3803/results/breakpoint_calling/sample_SA1256PP/destruct_breakpoints.csv.gz"


# #destruct lumpy
# SVc.consensus_destruct_lumpy(destruct, lumpy, "/juno/work/shah/abramsd/CODE/SV_bench/lumpy_destruct-consensus.csv", "/juno/work/shah/abramsd/CODE/SV_bench/lumpy_destruct-maches.csv", 500)
# #destruct svaba
# SVc.consensus_destruct_svaba(destruct, svaba, "/juno/work/shah/abramsd/CODE/SV_bench/destruct_svaba-consensus.csv", "/juno/work/shah/abramsd/CODE/SV_bench/destruct_svaba-matches.csv", 500)
# #lumpy svaba
# SVc.consensus_lumpy_svaba(lumpy, svaba, "/juno/work/shah/abramsd/CODE/SV_bench/lumpy_svaba-consensus.csv", "/juno/work/shah/abramsd/CODE/SV_bench/lumpy_svaba-.csv", 500)



consensus = pd.read_csv("consensus_matches/cleaned_match_lumpy_svaba.csv", sep="\t")

svaba = SVc.read_svaba(consensus)
not_matches = svaba[~svana.id.isin(consensus.id)]

print(svaba, not_matches)