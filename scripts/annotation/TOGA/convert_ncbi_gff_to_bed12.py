#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

import pandas as pd

from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output table with common name")
parser.add_argument("-a", "--allowed_types", action="store", dest="allowed_types", type=lambda s: s.split(","),
                    default=["gene", "pseudogene", "mRNA", "exon", "CDS"],
                    help="Comma-separated list of allowed annotations. Other types will be ignored. "
                         "Default: gene,pseudogene,mRNA,exon,CDS .")
args = parser.parse_args()

preprocessed_annotations = "{0}.preprocessed.bed".format(args.output_prefix)
print("Preprocessing started...")
with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(preprocessed_annotations, "w") as out_fd:
    out_fd.write("#scaffold\tstart\tend\tstrand\ttype\tid\tparent_id\tname\n")
    for line in in_fd:
        if line[0] == "#":
            continue

        line_list = line.strip().split("\t")
        #if line_list[2] not in args.allowed_types:
        #    continue
        description_dict = {key: value for key, value in map(lambda s: s.split("="), line_list[8].split(";"))}
        out_fd.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(line_list[0],
                                                                       int(line_list[3]) - 1,
                                                                       line_list[4],
                                                                       line_list[6],
                                                                       line_list[2],
                                                                       description_dict["ID"],
                                                                       description_dict["Parent"] if "Parent" in description_dict else ".",
                                                                       description_dict["Name"] if "Name" in description_dict else "."))

print("Preprocessing finished.")
print("Extracting isoforms...")
annotation_df = pd.read_csv(preprocessed_annotations, sep="\t", header=0, na_values=["."])
annotation_df[annotation_df["type"] == "mRNA"][["parent_id", "id"]].to_csv("{0}.isoforms.tab".format(args.output_prefix),
                                                                           sep="\t", header=False, index=False)
print("Extraction finished.")

parent_start_df = annotation_df[["id", "start"]]
parent_start_df.columns = pd.Index(["id", "parent_start"])

annotation_df = annotation_df.merge(parent_start_df, how='left', left_on="parent_id",
                                    right_on="id", suffixes=(None, "_parent"))
annotation_df["parent_start"] = annotation_df["parent_start"].astype("Int64")
annotation_df["parent_shift"] = annotation_df["start"] - annotation_df["parent_start"]

#[annotation_df.loc[parent_id, "start"] if parent_id != "." else 0 for parent_id in annotation_df["parent_id"]]
#annotation_df["parent_end"] = [annotation_df.loc[parent_id, "end"] if parent_id != "." else 0 for parent_id in annotation_df["parent_id"]]
#print(annotation_df)
annotation_df.to_csv("{0}.annotation_df.tab".format(args.output_prefix), sep="\t", header=True, index=True)


def exon_processing(df):
    #print("AAAAAA")
    #print(df["parent_shift"])
    #print(list(df["parent_shift"]))
    #print(",".join(map(str, list(df["parent_shift"]))) + ",")
    data_df = df.sort_values(by=["start", "end"])
    return pd.DataFrame.from_records([[len(data_df),
                                       ",".join(map(str, data_df["end"] - data_df["start"])) + ",",
                                       ",".join(map(str, data_df["parent_shift"])) + ","
                                       ]],
                                     columns=["exon_number", "exon_len_list", "exon_start_list"])
    #df["exon_number"] = df["exon_number"].astype("Int64")

    #return df


def cds_processing(df):
    data_df = df.sort_values(by=["start", "end"])
    return pd.DataFrame.from_records([[data_df["start"].iloc[0],
                                       data_df["end"].iloc[-1],
                                       ]],
                                     columns=["cds_start", "cds_end"])


print("Converting exon and CDS info...")

cds_df = annotation_df[annotation_df["type"] == "CDS"].groupby(["parent_id"]).apply(cds_processing)
exon_df = annotation_df[annotation_df["type"] == "exon"].groupby(["parent_id"]).apply(exon_processing)
cds_df["cds_start"] = cds_df["cds_start"].astype("Int64")
cds_df["cds_end"] = cds_df["cds_end"].astype("Int64")
exon_df["exon_number"] = exon_df["exon_number"].astype("Int64")

elements_df = pd.concat([cds_df, exon_df], axis=1)
elements_df.to_csv("{0}.exon.tab".format(args.output_prefix), sep="\t", header=True, index=True)

print("Conversion finished.")
#print(elements_df)

print("Preparing final files...")
merged_df = annotation_df.merge(elements_df, left_on="id", right_on="parent_id", how="right")
merged_df["score"] = 0
merged_df["score2"] = 0

bed12_columns = ["scaffold", "start", "end",
                 "transcript_id", "score", "strand",
                 "cds_start", "cds_end", "score2",
                 "exon_number", "exon_len_list", "exon_start_list"]
mRNA_index = merged_df["type"] == "mRNA"

merged_df = merged_df[["#scaffold", "start", "end",
                       "id", "score", "strand",
                       "cds_start", "cds_end", "score2",
                       "exon_number", "exon_len_list", "exon_start_list"]]

merged_df.columns = pd.Index(bed12_columns)

merged_df.to_csv("{0}.final.bed".format(args.output_prefix), sep="\t", index=False, header=False)
merged_df[mRNA_index].to_csv("{0}.final.mRNA.bed".format(args.output_prefix), sep="\t", index=False, header=False)
print("Finished...")
