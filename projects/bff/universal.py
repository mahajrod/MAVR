#!/usr/bin/env python
import sys
import argparse
from collections import OrderedDict
import numpy as np
import pandas as pd
from RouToolPa.Parsers.LAST import CollectionLast
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.Sequence import CollectionSequence

from RouToolPa.Routines import SequenceRoutines
from RouToolPa.Collections.General import IdList

parser = argparse.ArgumentParser()
"""
parser.add_argument("-g", "--genes", action="store", dest="genes", type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of names of genes to transfer.")

args = parser.parse_args()
# /home/mahajrod/tmp/annotation
gene_ids = args.genes
"""

parser.add_argument("-g", "--genes", action="store", dest="genes",  default=sys.stdin,
                    help="File with names of genes to transfer.")
parser.add_argument("-t", "--target_fasta", action="store", dest="target_fasta", required=True,
                    help="Fasta with target seqs.")
parser.add_argument("-q", "--query_fasta", action="store", dest="query_fasta", required=True,
                    help="Fasta with query seqs")
parser.add_argument("-l", "--last_tab", action="store", dest="last_tab", required=True,
                    help="LAST tab file")
parser.add_argument("-f", "--target_gff", action="store", dest="target_gff", required=True,
                    help="Target GFF")
parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", required=True,
                    help="Output directory")
args = parser.parse_args()

gene_ids = IdList(filename=args.genes)
print gene_ids

target_seq = CollectionSequence(
    in_file=args.target_fasta)
query_seq = CollectionSequence(
    in_file=args.query_fasta)

gff_coll = CollectionGFF(in_file=args.target_gff, # /home/skliver/df_spermatogenesis/homo_sapiens.gff
                         format="gff", parsing_mode="complete",
                         featuretype_separation=True,
                         scaffold_syn_dict=None)

last_coll = CollectionLast(in_file=args.last_tab, # /home/skliver/df_spermatogenesis/mustela_putorius_furo.to.homo_sapiens.R11.tab.gz
                           parsing_mode="complete")
output_dir = args.out_dir
SequenceRoutines.safe_mkdir(output_dir)

last_coll.records.index.name = "row"
cds_coords_df = gff_coll.records["CDS"][gff_coll.records["CDS"]["source"] == "BestRefSeq"][
    ["start", "end", "strand", "ID", "Parent"]]
cds_coords_df.index = cds_coords_df.index.droplevel(1)

# convert length column to end column. end and start coordinates are in python notation

last_coll.records["target_hit_len"] += last_coll.records["target_start"]
last_coll.records["query_hit_len"] += last_coll.records["query_start"]

last_coll.records.rename(columns={"target_hit_len": "target_end", "query_hit_len": "query_end"}, inplace=True)

# sort and index alignments

last_coll.records.sort_values(by=['target_id', "target_start", "query_id", "query_start"], inplace=True)

last_coll.records.set_index(["target_id", "query_id"], inplace=True, append="True", drop=False)


# last_coll.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
#                                                       names=("scaffold", "row"))

def parse_alignment_string(alignment_string, return_cum_sum=True):
    str_list = map(lambda s: s.split(":"), alignment_string.split(","))

    target_list = []
    query_list = []

    for entry in str_list:
        target_list.append(int(entry[0]))
        query_list.append(int(entry[0]) if len(entry) == 1 else int(entry[1]))
    target_list = np.array(target_list)
    query_list = np.array(query_list)

    if return_cum_sum:
        return target_list, query_list, np.cumsum(target_list), np.cumsum(query_list)
    else:
        return target_list, query_list


last_coll.records[["target_aln", "query_aln", "target_aln_cs", "query_aln_cs"]] = pd.DataFrame(
    map(parse_alignment_string, last_coll.records["alignment"]),
    columns=["target_aln", "query_aln", "target_aln_cs", "query_aln_cs"],
    index=last_coll.records.index)

source = "target"
if source == "target":
    source = "target"
elif source == "query":
    source = "query"
else:
    raise ValueError("ERROR!!! Unrecognized source! Only 'target' or 'query' are allowed.")

search = "query" if source == "target" else "target"

source_id = source + "_id"
search_id = search + "_id"

source_start = source + "_start"
search_start = search + "_start"

source_end = source + "_end"
search_end = search + "_end"

source_strand = source + "_strand"
search_strand = search + "_strand"

source_id = source + "_id"

source_aln_cs = source + "_aln_cs"
search_aln_cs = search + "_aln_cs"

start_aln_index = "start_aln_index"
start_source_position = "source_start_pos"
start_search_position = "transferred_start_pos"

end_aln_index = "end_aln_index"
end_source_position = "source_end_pos"
end_search_position = "transferred_end_pos"

start_source_position = "source_start_pos"
start_search_position = "transferred_start_pos"

source_seq_records, search_seq_records = (target_seq.records, query_seq.records) \
    if source == "target" else (query_seq.records, target_seq.records)


def transfer_start_func(s):
    if s[start_aln_index] == 0:

        return s[search_start] - s["Astart-start"]
    else:
        return s[search_start] + s[search_aln_cs][s[start_aln_index] - 1] - s["Astart-start"] - s[source_aln_cs][
            s[start_aln_index] - 1]


def transfer_end_func(s):
    if s[end_aln_index] == 0:
        return s[search_start] - s["Astart-end"]
    else:
        return s[search_start] + s[search_aln_cs][s[end_aln_index] - 1] - s["Astart-end"] - s[source_aln_cs][
            s[end_aln_index] - 1]


def transfer_coordinates(coordinates_tuple, verbose=False):
    def extact_source_seq():
        if coordinates_tuple[3] == "+":
            return source_seq_records[coordinates_tuple[0]][coordinates_tuple[1]: coordinates_tuple[2]]
        else:
            return SequenceRoutines.reverse_complement(
                source_seq_records[coordinates_tuple[0]][coordinates_tuple[1]: coordinates_tuple[2]])

    def extact_search_seq(row):
        # print row
        if row[source_strand] == "+":
            if row[search_strand] == "+":
                if coordinates_tuple[3] == "+":
                    return (search_seq_records[row[search_id]][row[start_search_position]:row[end_search_position]],
                            search_seq_records[row[search_id]][
                            row[start_search_position] - 2:row[start_search_position]],
                            search_seq_records[row[search_id]][row[end_search_position]:row[end_search_position] + 2],
                            row[start_search_position], row[end_search_position], "+")


                else:
                    return (SequenceRoutines.reverse_complement(
                        search_seq_records[row[search_id]][row[start_search_position]:row[end_search_position]]),
                            SequenceRoutines.reverse_complement(search_seq_records[row[search_id]][
                                                                row[end_search_position]:row[end_search_position] + 2]),
                            SequenceRoutines.reverse_complement(search_seq_records[row[search_id]][
                                                                row[start_search_position] - 2:row[
                                                                    start_search_position]]),
                            row[start_search_position], row[end_search_position], "-")
            else:

                scaffold_length = int(last_coll.query_scaffold_lengths["length"].loc[row[search_id]])
                new_start = scaffold_length - row[end_search_position]
                new_end = scaffold_length - row[start_search_position]

                if coordinates_tuple[3] == "+":
                    return (SequenceRoutines.reverse_complement(search_seq_records[row[search_id]][new_start:new_end]),
                            SequenceRoutines.reverse_complement(
                                search_seq_records[row[search_id]][new_end:new_end + 2]),
                            SequenceRoutines.reverse_complement(
                                search_seq_records[row[search_id]][new_start - 2:new_start]),
                            new_start, new_end, "-"
                            )

                else:
                    return (search_seq_records[row[search_id]][new_start:new_end],
                            search_seq_records[row[search_id]][new_start - 2:new_start],
                            search_seq_records[row[search_id]][new_end:new_end + 2],
                            new_start, new_end, "+")
        else:
            raise ValueError("ERROR!!! transfer from query not implemented yet")
        """    
        if ((coordinates_tuple[3] == "+") and (row[source_strand] == "+")) or ((coordinates_tuple[3] == "-") and (row[source_strand] == "-")):
            return search_seq_records[row[search_id]][row[start_search_position]:row[end_search_position]]
        else:
            return SequenceRoutines.reverse_complement(search_seq_records[row[search_id]][row[start_search_position]:row[end_search_position]])
        """

    last_col_chr = last_coll.records.xs(coordinates_tuple[0], level=source_id)

    coordinates_df = pd.DataFrame()
    coordinates_df[[source_start, source_end]] = last_col_chr[[source_start, source_end]]
    coordinates_df["Astart-start"] = last_col_chr[source_start] - coordinates_tuple[1]
    coordinates_df["Aend-start"] = last_col_chr[source_end] - coordinates_tuple[1]

    coordinates_df["Astart-end"] = last_col_chr[source_start] - coordinates_tuple[2]
    coordinates_df["Aend-end"] = last_col_chr[source_end] - coordinates_tuple[2]

    start_df = last_col_chr[(coordinates_df["Astart-start"] <= 0) & (coordinates_df["Aend-start"] > 0)].copy()
    if not start_df.empty:
        start_df["seq_id"] = coordinates_tuple[4]
        start_df["parent_id"] = coordinates_tuple[5]
        start_df["seq_strand"] = coordinates_tuple[3]
        start_df[["Astart-start", "Aend-start"]] = coordinates_df[["Astart-start", "Aend-start"]][
            (coordinates_df["Astart-start"] <= 0) & (coordinates_df["Aend-start"] > 0)]

        start_df[start_aln_index] = start_df.apply(lambda s: np.argmax(s[source_aln_cs] > -s["Astart-start"]), axis=1)
        start_df[start_source_position] = coordinates_tuple[1]
        start_df[start_search_position] = start_df.apply(transfer_start_func, axis=1)
        # start_df["start_index_offset"] = coordinates_df[[source_aln_cs]].apply(lambda s: np.argmax(s[source_aln_cs]>-s["Astart-start"]), axis=1)

    end_df = last_col_chr[(coordinates_df["Astart-end"] < 0) & (coordinates_df["Aend-end"] >= 0)].copy()

    if not end_df.empty:
        end_df["seq_id"] = coordinates_tuple[4]
        end_df["parent_id"] = coordinates_tuple[5]
        end_df["seq_strand"] = coordinates_tuple[3]
        end_df[["Astart-end", "Aend-end"]] = coordinates_df[["Astart-end", "Aend-end"]][
            (coordinates_df["Astart-end"] < 0) & (coordinates_df["Aend-end"] >= 0)]
        end_df[end_aln_index] = end_df.apply(lambda s: np.argmax(s[source_aln_cs] > -s["Astart-end"]), axis=1)
        end_df[end_source_position] = coordinates_tuple[2]
        end_df[end_search_position] = end_df.apply(transfer_end_func, axis=1)

    # print "START"
    # print start_df
    # print("\n\n")
    # print "END"
    # print end_df

    if (not end_df.empty) and (not start_df.empty):
        merged_df = pd.merge(start_df[["target_id", "target_strand", "query_id", "query_strand",
                                       "seq_id", "parent_id", "seq_strand",
                                       "source_start_pos", "transferred_start_pos"]],
                             end_df[[end_source_position,
                                     end_search_position]],
                             left_index=True, right_index=True)
    else:
        merged_df = pd.DataFrame(columns=[u'target_id', u'target_strand', u'query_id', u'query_strand', u'seq_id',
                                          u'parent_id', u'seq_strand', u'source_start_pos',
                                          u'transferred_start_pos', u'source_end_pos', u'transferred_end_pos',
                                          u'source_length', u'search_length', u'source_seq', u'search_seq',
                                          u'head_splice', u'trailing_splice', u'gff_start', u'gff_end',
                                          u'gff_strand'])

    if not merged_df.empty:
        # if seq_id not in merged_df.columns:
        #    merged_df["seq_id"] = coordinates_tuple[4]
        merged_df["source_length"] = merged_df[end_source_position] - merged_df[start_source_position]
        merged_df["search_length"] = merged_df[end_search_position] - merged_df[start_search_position]

        merged_df["source_seq"] = extact_source_seq()
        merged_df[
            ["search_seq", "head_splice", "trailing_splice", "gff_start", "gff_end", "gff_strand"]] = merged_df.apply(
            extact_search_seq, axis=1, result_type="expand")
    else:
        merged_df = pd.DataFrame([[coordinates_tuple[0], "",
                                               "", "",
                                               coordinates_tuple[4],
                                               coordinates_tuple[5], coordinates_tuple[3],
                                               coordinates_tuple[1], 0,
                                               coordinates_tuple[2], 0,
                                               0, 0,
                                               "", "", "", "",
                                               0, 0, ""]],
                                             columns=[u'target_id', u'target_strand',
                                                      u'query_id', u'query_strand',
                                                      u'seq_id',
                                                      u'parent_id', u'seq_strand',
                                                      u'source_start_pos', u'transferred_start_pos',
                                                      u'source_end_pos', u'transferred_end_pos',
                                                      u'source_length', u'search_length',
                                                      u'source_seq', u'search_seq',
                                                      u'head_splice', u'trailing_splice',
                                                      u'gff_start', u'gff_end', u'gff_strand'])
        merged_df["source_length"] = merged_df[end_source_position] - merged_df[start_source_position]
        merged_df["source_seq"] = extact_source_seq()
    if verbose:
        print("\n\n")
        print(merged_df)
    return start_df, end_df, merged_df


# gene_name-> gene_id -> mRNA id list
ids_dict = OrderedDict()


# genes_indel_ids = IdList(filename="/home/mahajrod/tmp/annotation/deletiom.insertion.indel.ids")
#all_ids = IdList(filename="/home/skliver/sp_spermatogenesis/all.ids")
for gene_name in gene_ids:
    ids_dict[gene_name] = OrderedDict()
    for gene_id in gff_coll.records["gene"][gff_coll.records["gene"]["Name"] == gene_name]["ID"]:
        ids_dict[gene_name][gene_id] = OrderedDict()
        for mrna_id in gff_coll.records["mRNA"][gff_coll.records["mRNA"]["Parent"] == gene_id]["ID"]:
            ids_dict[gene_name][gene_id][mrna_id] = OrderedDict({"CDS": [], "exon": []})

# gene_indel_gff_gene_ids = gff_coll.records["gene"][gff_coll.records["gene"]["Name"].isin(genes_indel_ids)]["ID"]
# gff_coll.records["mRNA"][gff_coll.records["mRNA"]["Parent"].isin(gene_indel_gff_ids)]["ID"]


SequenceRoutines.safe_mkdir(output_dir)
SequenceRoutines.recursive_mkdir(ids_dict, output_dir)

coords_df_dict = OrderedDict({"CDS": gff_coll.records["CDS"][["start", "end", "strand", "ID", "Parent"]],
                              "exon": gff_coll.records["exon"][["start", "end", "strand", "ID", "Parent"]]})
coords_df_dict["CDS"]["ID"] += "."
coords_df_dict["CDS"]["ID"] += map(str, np.arange(1, len(coords_df_dict["CDS"]) + 1))
for feature in coords_df_dict:
    coords_df_dict[feature].index = coords_df_dict[feature].index.droplevel(1)

for gene_name in ids_dict:
    # if gene_name in ("ADGRG2", "AR"):
    #    continue
    level = 0
    print(gene_name)
    for gene_id in ids_dict[gene_name]:
        level = 1
        print("\t" * level + gene_id)
        for mrna_id in ids_dict[gene_name][gene_id]:
            level = 2
            print("\t" * level + mrna_id)
            for feature in ["exon", "CDS"]:
                level = 3
                # print("\t"*level + feature)
                f_out_dir = "%s/%s/%s/%s/%s/" % (output_dir, gene_name, gene_id, mrna_id, feature)
                coords_df = coords_df_dict[feature][coords_df_dict[feature]["Parent"] == mrna_id]
                transfered_df_list = []
                for row_tuple in coords_df.itertuples(index=True):
                    level = 4
                    # print gff_coll.records[feature][gff_coll.records[feature]["Parent"] == mrna_id]
                    # print row_tuple
                    # print("\t"*level + "gff strand " +  row_tuple[3])

                    # print("\t"*level + row_tuple[4])
                    # print row_tuple
                    level = 5
                    start_df, end_df, merged_df = transfer_coordinates(row_tuple)

                    transfered_df_list.append(merged_df)
                    for table, extension in zip([start_df, end_df, merged_df], ["start", "end", "both"]):
                        """
                        if extension == "both":
                            #print table
                            print("\t"*level + "%s\t%i" % (extension, len(table) if int(table["search_length"].iloc[0]) !=0 else 0))
                        else:
                            print("\t"*level + "%s\t%i" % (extension, len(table) if not table.empty else 0))
                        """
                        if (not table.empty) and (extension != "both"):
                            # print extension
                            # print table.columns
                            table.to_csv("%s/%s.%s" % (f_out_dir, table["seq_id"][0], extension), sep="\t")
                else:
                    pd.concat(transfered_df_list).to_csv("%s/%s.%s.%s" % (f_out_dir, mrna_id, feature, "both"),
                                                         sep="\t")