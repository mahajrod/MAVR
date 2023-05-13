#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import numpy as np
import pandas as pd
from distinctipy import distinctipy


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--pcr_filelist", action="store", dest="pcr_filelist", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files containing pcr results.")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of labels for files containing pcr results."
                         "Must follow the same order as files.")
parser.add_argument("-d", "--loci_description_file", action="store", dest="loci_description_file", required=True,
                    help="File with description of amplified loci. Should look like:"
                         """primer_pair	source	min_length	max_length
Mf1.1	Basto et al, 2010	174	182
Mf1.11	Basto et al, 2010	219	234
Mf1.18	Basto et al, 2010	158	174
Mf1.3	Basto et al, 2010	221	237""")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")


args = parser.parse_args()

amplicon_info_df = pd.read_csv(args.loci_description_file, sep="\t", index_col="primer_pair")

amplicon_filedict = {label: file for label, file in zip(args.label_list, args.pcr_filelist)}

columns_to_use = ("amplicon_len", "HitName", "FP_ID", "FP_mismatches", "RP_ID",
                  "RP_mismatches", "Start", "End", "AmpliconSeq")

df_dict = {species: pd.read_csv(amplicon_filedict[species],
                                sep="\t", index_col=False, usecols=columns_to_use) for species in amplicon_filedict}
# ---- Init colors by source of loci ----
color_number = len(amplicon_info_df["source"])
colors = distinctipy.get_colors(color_number)
color_list = list(map(rgb_tuple_to_hex, colors))
# ----

for entry in df_dict:
    df_dict[entry].index.name = "row"
    df_dict[entry].set_index(keys=df_dict[entry]["FP_ID"].apply(lambda x: x.split("|")[0]), drop=False, inplace=True,
                             append=True)
    df_dict[entry].index.set_names(("row", "primer_pair"), inplace=True)
    df_dict[entry].index = df_dict[entry].index.swaplevel(0, 1)
    # df_dict[entry] = df_dict[entry][df_dict[entry]["amplicon_len"] <= 500]
    df_dict[entry]["max_mismatches"] = df_dict[entry][["FP_mismatches", "RP_mismatches"]].max(axis=1)
    df_dict[entry]["total_mismatches"] = df_dict[entry][["FP_mismatches", "RP_mismatches"]].sum(axis=1)

    df_dict[entry].sort_values(by=["primer_pair", "max_mismatches", "total_mismatches", "amplicon_len"], inplace=True)
    df_dict[entry][["max_mis_min_dist", "tot_mis_min_dist"]] = - df_dict[entry][
        ["max_mismatches", "total_mismatches"]].diff(periods=-1)

    df_dict[entry].loc[
        df_dict[entry].groupby('primer_pair').tail(1).index, ["max_mis_min_dist", "tot_mis_min_dist"]] = np.nan

    df_dict[entry]["min_len"] = 0
    df_dict[entry]["max_len"] = 0

    for pair in amplicon_info_df.index:
        if pair in df_dict[entry].index:
            df_dict[entry].loc[pair, "max_len"] = amplicon_info_df.loc[pair, "max_length"]
            df_dict[entry].loc[pair, "min_len"] = amplicon_info_df.loc[pair, "min_length"]

    df_dict[entry]["len_in_expected_interval"] = (0.5 * df_dict[entry]["min_len"] <= df_dict[entry]["amplicon_len"]) & (
                df_dict[entry]["amplicon_len"] <= 1.5 * df_dict[entry]["max_len"])

    df_dict[entry]["source"] = ""
    for primer_pair in amplicon_info_df.index:
        if primer_pair in df_dict[entry].index:
            df_dict[entry].loc[primer_pair, "source"] = amplicon_info_df.loc[primer_pair, "source"]

    df_dict[entry]["color"] = ""

    for source, color in zip(amplicon_info_df["source"].unique(), color_list):
                             #["red", "green", "blue", "orange", "violet"]):
        # print(list(amplicon_info_df[amplicon_info_df["Source"] == source].index))
        # print(df_dict[entry].index.unique(level=0))
        df_dict[entry].loc[(amplicon_info_df[amplicon_info_df["source"] == source].index & df_dict[entry].index.unique(
            level=0)), 'color'] = color

    amplicon_info_df[entry] = "NA"
    #print(df_dict[entry].index.get_level_values(level=0).unique())
    amplicon_info_df.loc[df_dict[entry].index.get_level_values(level=0).unique(), entry] = "A"

for species in amplicon_filedict:
    top_hits = df_dict[species].loc[df_dict[species].groupby('primer_pair').head(1).index]

    usable_STR_index = (top_hits["FP_ID"] != top_hits["RP_ID"]) & (top_hits["max_mismatches"] <= 2) & (
                top_hits["max_mis_min_dist"].isna() | (top_hits["max_mis_min_dist"] >= 2))
    possibly_usable_STR_index = (top_hits["FP_ID"] != top_hits["RP_ID"]) & (top_hits["max_mismatches"] <= 3) & (
                top_hits["max_mis_min_dist"].isna() | (top_hits["max_mis_min_dist"] >= 2) | (
                    top_hits["tot_mis_min_dist"] >= 3))
    possibly_usable_STR_index = (~usable_STR_index) & possibly_usable_STR_index

    all_usable_STR_index = usable_STR_index | possibly_usable_STR_index

    top_hits[usable_STR_index].reset_index(level=0).to_csv("{}.{}.usable.tsv".format(args.output_prefix, species),
                                                           sep="\t", header=True,
                                                           index=False)
    top_hits_usable_bed = top_hits[usable_STR_index].reset_index(level=0)[["HitName", "Start", "End", "primer_pair", "source", "color"]]
    top_hits_usable_bed.columns = pd.Index(["scaffold", "start", "end", "primer_pair", "source", "color"])
    top_hits_usable_bed.to_csv("{}.{}.usable.bed".format(args.output_prefix, species),
                               sep="\t", header=True, index=False)

    top_hits[possibly_usable_STR_index].reset_index(level=0).to_csv("{}.{}.possibly_usable.tsv".format(args.output_prefix, species), sep="\t",
                                                                    header=True, index=False)

    top_hits_possibly_usable_bed = top_hits[possibly_usable_STR_index].reset_index(level=0)[["HitName", "Start", "End", "primer_pair", "source", "color"]]
    top_hits_possibly_usable_bed.columns = pd.Index(["scaffold", "start", "end", "primer_pair", "source", "color"])
    top_hits_possibly_usable_bed.to_csv("{}.{}.possibly_usable.bed".format(args.output_prefix, species),
                                        sep="\t", header=True, index=False)

    top_hits[all_usable_STR_index].reset_index(level=0).to_csv("{}.{}.all_usable.tsv".format(args.output_prefix, species), sep="\t", header=True,
                                                               index=False)
    top_hits_all_usable_bed = top_hits[all_usable_STR_index].reset_index(level=0)[["HitName", "Start", "End", "primer_pair", "source", "color"]]
    top_hits_all_usable_bed.columns = pd.Index(["scaffold", "start", "end", "primer_pair", "source", "color"])
    top_hits_all_usable_bed.to_csv("{}.{}.all_usable.bed".format(args.output_prefix, species),
                                   sep="\t", header=True, index=False)
    amplicon_info_df.loc[
        all_usable_STR_index[all_usable_STR_index].index.get_level_values(level=0).unique(), species] = "U"

# ---- Color configuration ----
# -------- Color codes ----
green_hex = "#00FF00"
light_green_hex = "#90EE90"
light_blue_hex = "#ADD8E6"
light_yellow_hex = "#F4EA56"
light_orange_hex = "#FFD580"
light_red_hex = "#FF7377"
# --------

# ----

amplicon_info_df.to_csv("{0}.summary.tsv".format(args.output_prefix), sep="\t", index=True, header=True)

writer = pd.ExcelWriter('{0}.summary.xlsx'.format(args.output_prefix), engine='xlsxwriter')
workbook = writer.book

# -------- Coverage color formats --------
status_color_format_dict = {
                             "NA": workbook.add_format({'bg_color': light_red_hex}),  # not amplified
                             "A": workbook.add_format({'bg_color': light_yellow_hex}),  # amplified but failed criteria
                             #"Ok": workbook.add_format({'bg_color': light_blue_hex}),
                             #"High": workbook.add_format({'bg_color': light_green_hex}),
                             "U": workbook.add_format({'bg_color': green_hex})}  # useable - passe criteria

amplicon_info_df.to_excel(writer, sheet_name="summary")

header_row = 0
data_first_row = 1
data_last_row = data_first_row + len(amplicon_info_df) - 1
loci_column = 0
source_column = 1
min_length_column = 2
max_length_column = 3
data_first_column = 4
data_last_column = data_first_column + len(args.label_list) - 1

writer.sheets['summary'].freeze_panes(header_row, data_first_column)

for status in status_color_format_dict:
    writer.sheets['summary'].conditional_format(data_first_row,
                                                data_first_column,
                                                data_last_row,
                                                data_last_column,
                                                {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': status,
                                                 'format': status_color_format_dict[status]})
writer.sheets['summary'].set_column(loci_column, loci_column, 20)
writer.sheets['summary'].set_column(source_column, source_column, 30)
writer.sheets['summary'].set_column(min_length_column, max_length_column, 10)
writer.sheets['summary'].set_column(data_first_column, data_last_column, 20)
writer.save()
