#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import pandas as pd
import matplotlib as mpl
from functools import partial

from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Routines import SequenceRoutines
from RouToolPa.Tools.BLAST import MakeBLASTDb


def count_coverage_by_hits(df, start_column_name, end_column_name, final_col_prefix):
    #print(df)
    if len(df) <= 1:
        #return df[[start_column_name,end_column_name]]
        #print((df[end_column_name] - df[start_column_name]).iloc[0])
        return pd.DataFrame([(df[end_column_name] - df[start_column_name]).iloc[0]],
                            columns=[final_col_prefix + "_total_hit_len"])
    #print(df)
    start_col_idx = 5
    end_col_idx = 6
    row_list = [list(row) for row in df.iloc[0:1, :].itertuples(index=False)]
    #print(row_list)
    for row in df.iloc[1:, :].itertuples(index=False):
        if row[start_col_idx] > row_list[-1][end_col_idx]:
            row_list.append(row)
        else:
            if row[end_col_idx] > row_list[-1][end_col_idx]:
                row_list[-1][end_col_idx ] = row[end_col_idx]

    #return pd.DataFrame.from_records(row_list, columns=["query","target","target_strand", "query_len", start_column_name, end_column_name],)[[start_column_name,end_column_name]]
    tmp = pd.DataFrame.from_records(row_list, columns=["query","target", "query_len", "target_len", "target_strand",  start_column_name, end_column_name],)
    #print((tmp[end_column_name] - tmp[start_column_name]).sum())
    #return (tmp[end_column_name] - tmp[start_column_name]).sum()
    return pd.DataFrame([(tmp[end_column_name] - tmp[start_column_name]).sum()],
                        columns=[final_col_prefix + "_total_hit_len"])


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input_fasta", required=True,
                    help="Input fasta file ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-n", "--min_total_hit_len", action="store", dest="min_total_hit_len",
                    type=float, default=0.9,
                    help="Minimal threshold for total length of hits")
parser.add_argument("-m", "--max_total_hit_len", action="store", dest="max_total_hit_len",
                    type=float, default=1.1,
                    help="Maximal threshold for total length of hits")
args = parser.parse_args()

MakeBLASTDb.make_nucleotide_db(args.input_fasta, args.output_prefix, None,
                               output_file=args.output_prefix)

blast_hit_file = '%s.blastn.hists' % args.output_prefix

blast_cmd = 'blastn -outfmt "6 qaccver saccver qlen slen pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore" '
blast_cmd += ' -query %s ' % args.input_fasta
blast_cmd += ' -strand both '
blast_cmd += ' -db %s ' % args.output_prefix
blast_cmd += ' -out %s ' % blast_hit_file

os.system(blast_cmd)

blast_fmt6_header_list = ['query', 'target', "query_len", "target_len", 'pident', 'length', 'mismatch', 'gapopen',
                          'query_start', 'query_end', 'target_start', 'target_end', "target_strand", 'evalue', 'bitscore']

blast_df = pd.read_csv(blast_hit_file, sep="\t", names=blast_fmt6_header_list)

blast_df.loc[blast_df["target_strand"] == "minus", "target_start" ], blast_df.loc[blast_df["target_strand"] == "minus", "target_end"] = blast_df.loc[blast_df["target_strand"] == "minus", "target_end" ], blast_df.loc[blast_df["target_strand"] == "minus", "target_start" ]
blast_df["target_start"] = blast_df["target_start"] - 1
blast_df["query_start"] = blast_df["query_start"] - 1

blast_df["query_fraction"] = (blast_df["query_end"] - blast_df["query_start"]) / blast_df["query_len"]
blast_df["target_fraction"] = (blast_df["target_end"] - blast_df["target_start"]) / blast_df["target_len"]

query_sorted_blast_df = blast_df.sort_values(by=["query","target", "query_start"]) #[blast_df["query"] != blast_df["target"]]
target_sorted_blast_df = blast_df.sort_values(by=["query", "target", "target_start"]) #[blast_df["query"] != blast_df["target"]]

query_hit_len_df = query_sorted_blast_df[["query","target", "query_len","target_len","target_strand", "query_start", "query_end"]].groupby(by=["query", "target", "query_len","target_len",
                                                                                                           "target_strand"]).apply(partial(count_coverage_by_hits,
                                                                                                                                              start_column_name="query_start",
                                                                                                                                              end_column_name="query_end",
                                                                                                                                              final_col_prefix="query"))

query_hit_len_df.index = query_hit_len_df.index.droplevel(level=5)
#print(query_hit_len_df)
target_hit_len_df = target_sorted_blast_df[["query","target", "query_len","target_len","target_strand", "target_start", "target_end"]].groupby(by=["query", "target", "query_len","target_len",
                                                                                                              "target_strand"]).apply(partial(count_coverage_by_hits,
                                                                                                                                              start_column_name="target_start",
                                                                                                                                              end_column_name="target_end",
                                                                                                                                              final_col_prefix="target"))
target_hit_len_df.index = target_hit_len_df.index.droplevel(level=5)

hit_len_df = pd.concat([query_hit_len_df, target_hit_len_df], axis=1).reset_index(level=[2, 3], drop=False)
hit_len_df.to_csv('%s.hit_len.tab' % args.output_prefix, header=True, index=True, sep="\t")
hit_len_df["query_fraction"] = hit_len_df["query_total_hit_len"] / hit_len_df["query_len"]

hit_len_df["target_fraction"] = hit_len_df["target_total_hit_len"] / hit_len_df["target_len"]


hit_len_filtered_df = hit_len_df[(args.max_total_hit_len > hit_len_df["query_fraction"]) &
                                 (hit_len_df["query_fraction"] > args.min_total_hit_len) &
                                 (hit_len_df["query_fraction"] < args.max_total_hit_len) &
                                 (hit_len_df["target_fraction"] > args.min_total_hit_len)]
hit_len_filtered_df.to_csv('%s.hit_len.filtered.tab' % args.output_prefix, header=True, index=True, sep="\t")

singleton_set = set()
cluster_set = set()
for query_id in hit_len_filtered_df.index.get_level_values(level=0).unique():
    target_id_list = list(hit_len_filtered_df.loc[query_id].index.get_level_values(level=0))
    target_id_unique_list = list(hit_len_filtered_df.loc[query_id].index.get_level_values(level=0).unique())
    if len(target_id_list) != len(target_id_unique_list):
        raise ValueError("Some sequences have hits to both plus and minus strands of target. Check it.")
    if len(target_id_list) == 1:
        singleton_set.add(tuple(target_id_list))
    else:
        cluster_set.add(tuple(target_id_list))

print(singleton_set)
print(cluster_set)
with open("%s.clusters" % args.output_prefix, "w") as out_fd:
    for cluster_settttttt in (singleton_set, cluster_set):
        for cluster in cluster_settttttt:
            out_fd.write(",".join(cluster) + "\n")

indexed_query_sorted_blast_df = query_sorted_blast_df.set_index(["query", "target"])

hit_len_filtered_df.reset_index("target_strand", inplace=True)
modification_df = []
for cluster in singleton_set:

    # no modifications for singletons
    # seq_id, shift, revcomp
    modification_df.append([cluster[0], 0, False])

#print(modification_df)
for cluster in cluster_set:
    #detect most frequent strand
    first_seq_id = cluster[0]
    print(cluster)
    plus_strand_count = sum(hit_len_filtered_df.loc[first_seq_id ,"target_strand"] == "plus")
    minus_strand_count = sum(hit_len_filtered_df.loc[first_seq_id ,"target_strand"] == "minus")

    if minus_strand_count >= plus_strand_count:
        # change first seq if there more sequences in opposite strand
        first_seq_id = hit_len_filtered_df.loc[first_seq_id][hit_len_filtered_df.loc[first_seq_id]["target_strand"] == "minus"].index[0]
        #print(plus_strand_count, minus_strand_count)
    print(first_seq_id)
    for seq_id in cluster:
        if seq_id == first_seq_id:
            # no modifications for first(reference) seq in cluster
            # seq_id, shift, revcomp
            modification_df.append([first_seq_id, 0, False])
            continue
        if hit_len_filtered_df.loc[first_seq_id].loc[seq_id, "target_strand"] == "plus":
            modification_df.append([seq_id,
                                    indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["target_start"][0] - indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["query_start"][0],
                                    False])
        else:
            pass
            modification_df.append([seq_id,
                                    indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["target_end"][0] + indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["query_start"][0],
                                    True])

        print("\t" + str(modification_df[-1]))


modification_df = pd.DataFrame.from_records(modification_df, columns=["seq_id", "shift", "rev_com"], index="seq_id")
modification_df.to_csv('%s.modifications.tab' % args.output_prefix, header=True, index=True, sep="\t")

fasta_col = CollectionSequence(args.input_fasta, format="fasta", parsing_mode="parse")
for seq_id in modification_df.index:
    shift = modification_df.loc[seq_id, "shift"]
    rev_com = modification_df.loc[seq_id, "rev_com"]
    #if (modification_df.loc[seq_id, "shift"] == 0) and (not modification_df.loc[seq_id, "revcom"]):
    #    continue
    if shift != 0:
        fasta_col.records[seq_id] = fasta_col.records[seq_id][shift:] + fasta_col.records[seq_id][:shift]
    if rev_com:
        fasta_col.records[seq_id] = SequenceRoutines.reverse_complement(fasta_col.records[seq_id])

fasta_col.write("%s.rotated.fasta" % args.output_prefix)




