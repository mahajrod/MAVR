#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

import pandas as pd

from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", required=True,
                    help="Output table with common name")
parser.add_argument("-a", "--allowed_types", action="store", dest="allowed_types", type=lambda s: s.split(","),
                    default=["gene", "pseudogene", "mRNA", "exon", "CDS"],
                    help="Comma-separated list of allowed annotations. Other types will be ignored. "
                         "Default: gene,pseudogene,mRNA,exon,CDS .")
args = parser.parse_args()

preprocessed_annotations = "{0}.preprocessed.bed".format(args.output_prefix)

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(preprocessed_annotations, "w") as out_fd:
    out_fd.write("#scaffold\tstart\tend\tstrand\ttype\tid\tparent_id\tname\n")
    for line in in_fd:
        if line[0] == "#":
            continue

        line_list = line.strip().split("\t")
        if line_list[2] not in args.allowed_types:
            continue
        description_dict = {key: value for key, value in map(lambda s: s.split("="), line_list[8].split(";"))}
        out_fd.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(line_list[0],
                                                                       int(line_list[3]) - 1,
                                                                       line_list[4],
                                                                       line_list[6],
                                                                       line_list[2],
                                                                       description_dict["ID"],
                                                                       description_dict["Parent"] if "Parent" in description_dict else ".",
                                                                       description_dict["Name"] if "Name" in description_dict else "."))


annotation_df = pd.read_csv(preprocessed_annotations, sep="\t", header=0)
annotation_df[annotation_df["type"] == "mRNA"][["parent_id", "id"]].to_csv("{0}.isoforms.tab".format(args.output_prefix),
                                                                           sep="\t", header=False, index=False)
annotation_df.set_index("id", inplace=True)
annotation_df["parent_start"] = [annotation_df.loc[parent_id, "start"] if parent_id != "." else 0 for parent_id in annotation_df["parent_id"]]
annotation_df["parent_end"] = [annotation_df.loc[parent_id, "end"] if parent_id != "." else 0 for parent_id in annotation_df["parent_id"]]
print(annotation_df)


def exon_processing(df):
    return pd.DataFrame.from_records([[len(df),
                                       ",".join(map(str, df["end"] - df["start"])) + ",",
                                       ",".join(map(str, df["start"] - df["parent_start"])) + ","
                                       ]],
                                     columns=["exon_number", "exon_len_list", "exon_start_list"])


def cds_processing(df):
    return pd.DataFrame.from_records([[df["start"][0],
                                       df["end"][0],
                                       ]],
                                     columns=["cds_start", "cds_end",])


cds_df = annotation_df[annotation_df["type"] == "CDS"].groupby(["parent_id"]).apply(cds_processing)
exon_df = annotation_df[annotation_df["type"] == "exon"].groupby(["parent_id"]).apply(exon_processing)

print(cds_df)
print(exon_df)

