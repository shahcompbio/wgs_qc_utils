from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree
from scgenome import csvutils 
from wgs_qc_utils.reader import read_variant_calls
from wgs_qc_utils.reader import read_variant_calls
from wgs.utils import vcfutils
import pandas as pd

def build_interval_tree(data):
    '''
    load lumpy confidence intervals into tree
    '''
    itree = defaultdict(IntervalTree)

    for val in data:
        chrom, pos, ci, brk_type, id = val

        ci = ci.split(',')

        start = pos + float(ci[0])
        end = pos + float(ci[1]) + 1

        itree[chrom].addi(start, end, {"strand":brk_type, "id":id})

    return itree


def check_olp(interval_tree, chrom, pos, destruct_strand, id, svaba=False):
    
    overlap = interval_tree[chrom][pos]
    is_match =False
    if svaba:
        if overlap:
            # print(overlap, chrom, pos)
            is_match=True

    overlapping_strands = [v[2]["strand"] for v in overlap]

    if destruct_strand in overlapping_strands:
        is_match=True


    return is_match


def get_matches(interval_tree, chrom, pos, destruct_strand, id, svaba=False):
   matches = interval_tree[chrom][pos]
   matches = [{"chrom":chrom, "pos":match[0]+500, "strand":match[2]["strand"], "id":match[2]["id"], "distance": abs(pos-(match[0]+500)) } for match in matches]
   get_closest = pd.DataFrame(matches) 
   return get_closest[get_closest.pos == get_closest.pos.min()].to_dict(orient="records")[0]

def load_data(infile):

    csv_input = csvutils.CsvInput(infile)
    return csv_input.read_csv()



def read_svaba(f):
    data = read_variant_calls.read_full_slow(f, gz=True)[["ID", "CHROM", "POS", "MATEID"]]
    data["mate_id"] =  data.MATEID.apply(lambda id: id.split(":")[0])
    data["mate_num"] =  data.MATEID.apply(lambda id: id.split(":")[1])
    data_ones = data[data.mate_num == "1"]
    data_twos = data[data.mate_num == "2"]
    data = data_twos.merge(data_ones, on="mate_id")
    data = data.rename(columns={"mate_id": "id", "CHROM_x": "chromosome_1", "CHROM_y": "chromosome_2",
        "POS_x": "position_1", "POS_y": "position_2", 
        "mate_num_x": "mate_num_1", "mate_num_y": "mate_num_2"})
    return data[["id", "chromosome_1", "position_1", "chromosome_2", "position_2", "mate_num_1", "mate_num_2"]]



def load_lumpy_into_tree(lumpy_df, confidence_interval=None, svaba=False):
    # Add lumpy breakpoint id to each zipped entry
    if svaba: 
        strand_1 = ["?"] * len(lumpy_df)
        strand_2 = strand_1
    else:
        strand_1 = lumpy_df.strand_1
        strand_2 = lumpy_df.strand_2
    if confidence_interval:
        confidence_interval = '-{},{}'.format(confidence_interval, confidence_interval)
        data = list(
            zip(lumpy_df.chromosome_1, lumpy_df.position_1, [confidence_interval] * len(lumpy_df), strand_1, lumpy_df.id))
        data += list(
            zip(lumpy_df.chromosome_2, lumpy_df.position_2, [confidence_interval] * len(lumpy_df), strand_2, lumpy_df.id))
    else:
        data = list(zip(lumpy_df.chromosome_1, lumpy_df.position_1, lumpy_df['CIPOS'], lumpy_df['strand_1'], lumpy_df.id))
        data += list(zip(lumpy_df.chromosome_2, lumpy_df.position_2, lumpy_df['CIEND'], lumpy_df['strand_2'], lumpy_df.id))

    intervaltree = build_interval_tree(data)

    return intervaltree


def filter_destruct_on_lumpy(destruct, lumpy_tree, svaba=False):

    if not svaba:
        destruct = destruct[
            destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_1'], x['position_1'], x['strand_1'], x["id"], svaba=svaba),
                        axis=1)]
        destruct = destruct[
            destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_2'], x['position_2'], x['strand_2'], x["id"], svaba=svaba),
                        axis=1)]
        print(destruct, "EARLY")
    if svaba:
        destruct = destruct[
            destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_1'], x['position_1'],"?",x["id"], svaba=svaba),
                        axis=1)]
        print(destruct, "EARLY")

        destruct = destruct[
            destruct.apply(lambda x: check_olp(lumpy_tree, x['chromosome_2'], x['position_2'], "?",x["id"], svaba=svaba),
                        axis=1)]      

    if destruct.empty:    
        return destruct, destruct, destruct
    
    else:
        matches = destruct
        if svaba:
            matches["matches_2"] = matches.apply(lambda x: 
                get_matches(lumpy_tree, x['chromosome_2'], x['position_2'], "?", x["id"], svaba=svaba), axis=1)
            matches["matches_1"] = matches.apply(lambda x: 
                get_matches(lumpy_tree, x['chromosome_1'], x['position_1'], "?", x["id"], svaba=svaba), axis=1)
        else:
            matches["matches_2"] = matches.apply(lambda x: 
                get_matches(lumpy_tree, x['chromosome_2'], x['position_2'], x["strand_2"], x["id"], svaba=svaba), axis=1)
            matches["matches_1"] = matches.apply(lambda x: 
                get_matches(lumpy_tree, x['chromosome_1'], x['position_1'], x["strand_1"], x["id"], svaba=svaba), axis=1)

        matches = matches[(matches.matches_1 !=set()) & (matches.matches_2!=set())]
        cleaned_matches=pd.DataFrame()
        cleaned_matches["id_1"]=  matches.id
        cleaned_matches["id_2_a"] = matches.matches_1.apply(lambda x : x["id"])
        cleaned_matches["id_2_a_distance"] = matches.matches_1.apply(lambda x : x["distance"])
        cleaned_matches["id_2_b"] = matches.matches_2.apply(lambda x : x["id"])
        cleaned_matches["id_2_b_distance"] = matches.matches_2.apply(lambda x : x["distance"])

        # cleaned_matches[["id_1", "id_1_a", "id_1_a_distance", "id_2_a", "id_2_a_distance"]] = matches.apply(lambda x: (x["id"], x["matches_1"]["id"], x["matches_1"]["distance"], x["matches_2"]["id"],x["matches_2"]["distance"]), axis=1)
        return destruct, matches, cleaned_matches
        
def write(data, outfile):
    data.to_csv(outfile, index=False)


def consensus_destruct_lumpy(destruct_infile, lumpy_infile, consensus, matches, confidence_interval=None):
    destruct = load_data(destruct_infile)
    destruct = destruct.astype({'chromosome_1': str, 'chromosome_2': str})
    destruct = destruct.rename(columns={"prediction_id": "id"})

    lumpy = load_data(lumpy_infile).rename(columns= {
        "chrom1": "chromosome_1", "chrom2": "chromosome_2", 
        "start1": "position_1", "start2": "position_2", 
        "strand1": "strand_1", "strand2": "strand_2", "breakpoint_id": "id"})

    lumpy.astype({'chromosome_1': str, 'chromosome_2': str}, inplace=True)

    lumpy = load_lumpy_into_tree(lumpy, confidence_interval=confidence_interval)
    print(len(lumpy))
    destruct,matches_output, cleaned_matches = filter_destruct_on_lumpy(destruct, lumpy)

    write(destruct, consensus)
    matches_output.to_csv(matches, sep="\t", index=False)
    cleaned_matches.to_csv("cleaned_match_destruct_lumpy.csv", sep="\t", index=False)


def consensus_destruct_svaba(destruct_infile, svaba_infile, consensus, matches, confidence_interval=None):
    destruct = load_data(destruct_infile)
    destruct = destruct.astype({'chromosome_1': str, 'chromosome_2': str})
    destruct = destruct.rename(columns={"prediction_id": "id"})
    svaba = read_svaba(svaba_infile)
    svaba = load_lumpy_into_tree(svaba, confidence_interval=confidence_interval, svaba=True)
    destruct, matches_output, cleaned_matches = filter_destruct_on_lumpy(destruct, svaba, svaba=True)

    write(destruct, consensus)


    matches_output.to_csv(matches, sep="\t", index=False)
    cleaned_matches.to_csv("cleaned_match_destruct_svaba.csv", sep="\t", index=False)

def consensus_lumpy_svaba(lumpy_infile, svaba_infile, consensus, matches, confidence_interval=None):
    lumpy = load_data(lumpy_infile).rename(columns= {
        "chrom1": "chromosome_1", "chrom2": "chromosome_2", 
        "start1": "position_1", "start2": "position_2", 
        "strand1": "strand_1", "strand2": "strand_2", "breakpoint_id":"id"})

    lumpy.astype({'chromosome_1': str, 'chromosome_2': str}, inplace=True)

    svaba = read_svaba(svaba_infile)
    print(svaba)
    lumpy = load_lumpy_into_tree(lumpy, confidence_interval=confidence_interval, svaba=True)
    svaba, matches_output, cleaned_matches = filter_destruct_on_lumpy(svaba, lumpy, svaba=True)

    write(svaba, consensus)
    matches_output.to_csv(matches, sep="\t", index=False)
    print(cleaned_matches)
    cleaned_matches.to_csv("cleaned_match_lumpy_svaba.csv", sep="\t", index=False)