#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from copy import deepcopy
from functools import partial
from collections.abc import Iterable
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Routines import SequenceRoutines
from RouToolPa.Tools.BLAST import MakeBLASTDb


def count_coverage_by_hits(df, start_column_name, end_column_name, final_col_prefix):
    #print("OOOO")
    #print(df)
    #print(type(df))
    tmp_df = df[[start_column_name, end_column_name]].reset_index(drop=True)
    #print(tmp_df)
    if len(tmp_df) <= 1:
        #print("AAAAAAAA")
        #print(df)
        #return df[[start_column_name,end_column_name]]
        #print((df[end_column_name] - df[start_column_name]).iloc[0])
        return pd.DataFrame([(tmp_df[end_column_name] - tmp_df[start_column_name]).iloc[0]],
                            columns=[final_col_prefix + "_total_hit_len"])
    #print("BBBBBBBB")

    #print(df)
    start_col_idx = 0 #5
    end_col_idx = 1 #6
    row_list = [list(row) for row in tmp_df.iloc[0:1, :].itertuples(index=False)]
    #print(row_list)
    #print(row_list)
    for row in tmp_df.iloc[1:, :].itertuples(index=False):
        tmp_row = list(row)
        if tmp_row[start_col_idx] > row_list[-1][end_col_idx]:
            row_list.append(tmp_row)
        else:
            if tmp_row [end_col_idx] > row_list[-1][end_col_idx]:
                #print("CCCCCC")
                #print(tmp_row )
                #print(row_list)
                row_list[-1][end_col_idx] = tmp_row[end_col_idx]

    #return pd.DataFrame.from_records(row_list, columns=["query","target","target_strand", "query_len", start_column_name, end_column_name],)[[start_column_name,end_column_name]]
    #tmp = pd.DataFrame.from_records(row_list, columns=["query","target", "query_len", "target_len", "target_strand",  start_column_name, end_column_name],)
    tmp = pd.DataFrame.from_records(row_list, columns=[start_column_name, end_column_name],)
    #print((tmp[end_column_name] - tmp[start_column_name]).sum())
    #return (tmp[end_column_name] - tmp[start_column_name]).sum()
    return pd.DataFrame([(tmp[end_column_name] - tmp[start_column_name]).sum()],
                        columns=[final_col_prefix + "_total_hit_len"])


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input_fasta", required=True,
                    help="Input fasta file ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-m", "--max_difference", action="store", dest="max_difference",
                    type=float, default=0.1,
                    help="Maximal threshold for difference between sequences "
                         "(fraction of sequence uncovered by blast hits). Default: 0.1")
parser.add_argument("-e", "--evalue", action="store", dest="evalue",
                    type=float, default=0.001,
                    help="Maximal threshold for e-value (blast). Default: 0.001")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.05,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")

parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float, default=0.98,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")

parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float, default=0.90,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float, default=0.05,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")

args = parser.parse_args()

MakeBLASTDb.make_nucleotide_db(args.input_fasta, args.output_prefix, None,
                               output_file=args.output_prefix)

blast_hit_file = '%s.blastn.hits' % args.output_prefix

blast_cmd = 'blastn -outfmt "6 qaccver saccver qlen slen pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore" '
blast_cmd += ' -evalue %f' % args.evalue
blast_cmd += ' -query %s ' % args.input_fasta
blast_cmd += ' -strand both '
blast_cmd += ' -db %s ' % args.output_prefix
blast_cmd += ' -out %s ' % blast_hit_file

os.system(blast_cmd)

blast_fmt6_header_list = ['query', 'target', "query_len", "target_len", 'pident', 'length', 'mismatch', 'gapopen',
                          'query_start', 'query_end', 'target_start', 'target_end', "target_strand", 'evalue', 'bitscore']

blast_df = pd.read_csv(blast_hit_file, sep="\t", names=blast_fmt6_header_list)
fasta_col = CollectionSequence(args.input_fasta, format="fasta", parsing_mode="parse")

blast_df.loc[blast_df["target_strand"] == "minus", "target_start"], blast_df.loc[blast_df["target_strand"] == "minus", "target_end"] = blast_df.loc[blast_df["target_strand"] == "minus", "target_end"], blast_df.loc[blast_df["target_strand"] == "minus", "target_start"]
blast_df["target_start"] = blast_df["target_start"] - 1
blast_df["query_start"] = blast_df["query_start"] - 1

blast_df["query_fraction"] = (blast_df["query_end"] - blast_df["query_start"]) / blast_df["query_len"]
blast_df["target_fraction"] = (blast_df["target_end"] - blast_df["target_start"]) / blast_df["target_len"]

query_sorted_blast_df = blast_df.sort_values(by=["query", "target", "query_start"]) #[blast_df["query"] != blast_df["target"]]
target_sorted_blast_df = blast_df.sort_values(by=["query", "target", "target_start"]) #[blast_df["query"] != blast_df["target"]]

#
#def test (df):
#    print(df)
#    return 0
#

query_sorted_blast_df[["query","target", "query_len","target_len","target_strand", "query_start", "query_end"]].to_csv("aaaa", sep="\t", index=True, header=True)
#query_sorted_blast_df[["query","target", "query_len","target_len","target_strand", "query_start", "query_end"]].groupby(by=["query", "target", "query_len","target_len",
#                                                                                                           "target_strand"]).apply(test)

query_hit_len_df = query_sorted_blast_df[["query", "target", "query_len", "target_len", "target_strand", "query_start", "query_end"]].groupby(by=["query", "target", "query_len", "target_len",
                                                                                                           "target_strand"]).apply(partial(count_coverage_by_hits,
                                                                                                                                              start_column_name="query_start",
                                                                                                                                              end_column_name="query_end",
                                                                                                                                              final_col_prefix="query"))

query_hit_len_df.index = query_hit_len_df.index.droplevel(level=5)
#print(query_hit_len_df)
target_hit_len_df = target_sorted_blast_df[["query", "target",
                                            "query_len", "target_len",
                                            "target_strand", "target_start",
                                            "target_end"]].groupby(by=["query", "target",
                                                                       "query_len", "target_len",
                                                                       "target_strand"]).apply(partial(count_coverage_by_hits,
                                                                       start_column_name="target_start",
                                                                       end_column_name="target_end",
                                                                       final_col_prefix="target"))
target_hit_len_df.index = target_hit_len_df.index.droplevel(level=5)

hit_len_df = pd.concat([query_hit_len_df, target_hit_len_df], axis=1).reset_index(level=[2, 3], drop=False)
hit_len_df.to_csv('%s.hit_len.tab' % args.output_prefix, header=True, index=True, sep="\t")
hit_len_df["query_fraction"] = hit_len_df["query_total_hit_len"] / hit_len_df["query_len"]

hit_len_df["target_fraction"] = hit_len_df["target_total_hit_len"] / hit_len_df["target_len"]

print("Selecting main strand for each pair of target and query...")

row_list = [list(row) for row in hit_len_df.reset_index(drop=False).iloc[0:1, :].itertuples(index=False)]
#print(hit_len_df)
for row in hit_len_df.iloc[1:, :].reset_index(drop=False).itertuples(index=False):
    if (row[0] == row_list[-1][0]) and (row[1] == row_list[-1][1]):
        if (row[-1] - row_list[-1][-1] + row[-2] - row_list[-1][-2]) > 0:
            row_list[-1] = deepcopy(row)
    else:
        row_list.append(row)

hit_len_df = pd.DataFrame.from_records(row_list, columns=["query", "target", "target_strand",
                                                          "query_len", "target_len",
                                                          "query_total_hit_len", "target_total_hit_len",
                                                          "query_fraction", "target_fraction"],
                                       index=["query", "target", "target_strand"])
hit_len_df.to_csv('%s.hit_len.strand_filtered.tab' % args.output_prefix, header=True, index=True, sep="\t")

seq_number = len(fasta_col.records)
seq_list = sorted(fasta_col.records.keys())
seq_index_dict = {seq_list[i]: i for i in range(0, seq_number)}

#distance_matrix = np.ones((seq_number, seq_number))
#hit_len_df_index = hit_len_df.index.droplevel([2, 3, 4])
#for query_id in fasta_col.records:
#    for target_id in fasta_col.records:
#        if (query_id, target_id) in hit_len_df_index:
#            distance_matrix[seq_index_dict[query_id]][seq_index_dict[target_id]] = max(np.abs(1 - hit_len_df.loc[query_id].loc[target_id]["query_total_hit_len"]),
#                                                                                       np.abs(1 - hit_len_df.loc[query_id].loc[target_id]["target_total_hit_len"]))
print("Clustering...")
distance_matrix = np.ones(seq_number * (seq_number - 1) // 2, dtype=float)
distance_table = np.ones((seq_number, seq_number), dtype=np.float64)
#print(distance_table.shape)
#print(type(distance_table))
for i in range(0, seq_number):
    distance_table[i][i] = 0
#print(hit_len_df.index)
hit_len_df_index = hit_len_df.index.droplevel(2)

for target_index in range(0, seq_number):  # j
    target_id = seq_list[target_index]
    for query_index in range(0, target_index):  # i
        query_id = seq_list[query_index]
        distance_matrix_index = seq_number * query_index + target_index - ((query_index + 2) * (query_index + 1)) // 2
        if ((query_id, target_id) in hit_len_df_index) and ((target_id, query_id) in hit_len_df_index):
            distance_table[query_index][target_index] = max(np.abs(1 - hit_len_df.loc[query_id].loc[target_id]["query_fraction"][0]),
                                                            np.abs(1 - hit_len_df.loc[query_id].loc[target_id]["target_fraction"][0]))
            distance_table[target_index][query_index] = max(np.abs(1 - hit_len_df.loc[target_id].loc[query_id]["query_fraction"][0]),
                                                            np.abs(1 - hit_len_df.loc[target_id].loc[query_id]["target_fraction"][0]))
            distance_matrix[distance_matrix_index] = max(distance_table[query_index][target_index],
                                                         distance_table[target_index][query_index])
#print(distance_table)
distance_df = pd.DataFrame(distance_table, index=seq_list, columns=seq_list)
distance_df.index.name = "seq_id"
distance_df.to_csv("%s.distance.tab" % args.output_prefix, sep="\t", index=True, header=True)
np.savetxt("%s.distance" % args.output_prefix, distance_matrix, fmt='%2.5f', delimiter="\t")
linkage = hierarchy.linkage(distance_matrix, method="complete")
fig = plt.figure(figsize=(25, 10))
dn = hierarchy.dendrogram(linkage, labels=seq_list)
plt.axhline(0.1, color="green", linestyle="dotted")
plt.axhline(0.2, color="blue", linestyle="dotted")
plt.axhline(0.3, color="orange", linestyle="dotted")
plt.axhline(0.5, color="red", linestyle="dotted")
plt.xlabel("Sequence")
plt.ylabel("Distance")
plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right, bottom=args.subplots_adjust_bottom,
                    top=args.subplots_adjust_top)

for ext in "png", "svg":
    plt.savefig("{0}.clustering.{1}".format(args.output_prefix, ext))
cluster_array = hierarchy.fcluster(linkage, args.max_difference, criterion='distance')

#print(cluster_array)
cluster_dict = {}
for seq_index in range(0, len(cluster_array)):
    cluster_index = cluster_array[seq_index]
    cluster_label = "cluster_%i" % cluster_index
    if cluster_label not in cluster_dict:
        cluster_dict[cluster_label] = [seq_list[seq_index]]
    else:
        cluster_dict[cluster_label].append(seq_list[seq_index])

print("Sequence indexes:")
for seq_id in seq_index_dict:
    print("\t{0}: {1}".format(seq_id, seq_index_dict[seq_id]))

#print(cluster_dict)

print("Clusters:")
for cluster_label in cluster_dict:
    print("\t{0}\t{1}".format(cluster_label, ",".join(cluster_dict[cluster_label])))


with open("%s.clusters" % args.output_prefix, "w") as out_fd:
    for cluster_label in cluster_dict:
        out_fd.write("{0}\t{1}\n".format(cluster_label, ",".join(cluster_dict[cluster_label])))


hit_len_df.reset_index("target_strand", inplace=True)

indexed_query_sorted_blast_df = query_sorted_blast_df.set_index(["query", "target"])

modification_df = []
for cluster_label in cluster_dict:
    if len(cluster_dict[cluster_label]) == 1:
        modification_df.append([cluster_dict[cluster_label][0], 0, False])
        continue
    #detect most frequent strand
    first_seq_id = cluster_dict[cluster_label][0]
    print("Cluster:\n\tid:\t%s" % str(cluster_label))
    print("\telements:\t%s" % str(",".join(cluster_dict[cluster_label])))
    plus_strand_count = sum(hit_len_df.loc[first_seq_id].loc[cluster_dict[cluster_label]]["target_strand"] == "plus")
    minus_strand_count = sum(hit_len_df.loc[first_seq_id].loc[cluster_dict[cluster_label]]["target_strand"] == "minus")
    print("\tplus strand: %i" % plus_strand_count)
    print("\tminus strand: %i" % minus_strand_count)
    print("\tInitial reference seq id:\t%s" % first_seq_id)
    if minus_strand_count >= plus_strand_count:
        # change first seq if there more sequences in opposite strand
        first_seq_id = hit_len_df.loc[first_seq_id].loc[cluster_dict[cluster_label]][hit_len_df.loc[first_seq_id].loc[cluster_dict[cluster_label]]["target_strand"] == "minus"].index[0]
        #print(plus_strand_count, minus_strand_count)
    print("\tNew reference seq id:\t%s" % first_seq_id)
    for seq_id in cluster_dict[cluster_label]:
        #print(seq_id)
        if seq_id == first_seq_id:
            # no modifications for first(reference) seq in cluster
            # seq_id, shift, revcomp
            shift = 0
            rev_com = False
            #modification_df.append([first_seq_id, shift, rev_com])
        else:
            #print(hit_len_df.loc[first_seq_id])
            #print(indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id])
            #print(indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["target_end"])
            #print(indexed_query_sorted_blast_df.loc[first_seq_id].loc[seq_id]["query_start"])
            if hit_len_df.loc[first_seq_id].loc[seq_id, "target_strand"] == "plus":
                shift = indexed_query_sorted_blast_df.loc[first_seq_id]["target_start"].loc[[seq_id]][0] - indexed_query_sorted_blast_df.loc[first_seq_id]["query_start"].loc[[seq_id]][0]
                if isinstance(shift, Iterable):
                    if (len(shift)) > 1:
                        raise ValueError("ERROR! strange issue with %s" % seq_id)
                    elif isinstance(shift, Iterable):
                        shift = shift[0]
                #print(indexed_query_sorted_blast_df.loc[first_seq_id]["target_start"].loc[[seq_id]])
                #print(indexed_query_sorted_blast_df.loc[first_seq_id]["query_start"].loc[[seq_id]])
                #print(shift)
                rev_com = False
            else:
                #print(indexed_query_sorted_blast_df.loc[first_seq_id]["target_end"].loc[[seq_id]])
                #print(indexed_query_sorted_blast_df.loc[first_seq_id]["query_start"].loc[[seq_id]])
                shift = indexed_query_sorted_blast_df.loc[first_seq_id]["target_end"].loc[[seq_id]][0] + indexed_query_sorted_blast_df.loc[first_seq_id]["query_start"].loc[[seq_id]][0],
                if isinstance(shift, Iterable):
                    if (len(shift)) > 1:
                        raise ValueError("ERROR! Strange issue with %s happened while calculating modifications" % seq_id)
                    elif isinstance(shift, Iterable):
                        shift = shift[0]
                #print(shift)
                rev_com = True
        modification_df.append([seq_id, shift, rev_com])
        #print(modification_df[-1])

        #print("\t%s\t%i\t%s" % (seq_id, shift, str(rev_com)))

        #print("\t" + str(modification_df[-1]))
    print("\n")

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
for cluster_label in cluster_dict:
    seq_file = "{0}.rotated.{1}.fasta".format(args.output_prefix, cluster_label)
    aln_file = "{0}.rotated.{1}.aln.fasta".format(args.output_prefix, cluster_label)
    fasta_col.write(seq_file,
                    whitelist=cluster_dict[cluster_label])

    mafft_cmd = "mafft --quiet --thread 4 --maxiterate 1000 --auto --anysymbol {0} > {1}".format(seq_file, aln_file)
    os.system(mafft_cmd)


