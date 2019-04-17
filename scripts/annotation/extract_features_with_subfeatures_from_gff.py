#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
# TODO: NOT FINISHED
import argparse
from RouToolPa.Collections.General import IdList

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output .gff file")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids of genes to extract")
parser.add_argument("-w", "--write_comments", action="store_true", dest="write_comments",
                    help="Write comments to output")

args = parser.parse_args()

feature_id_list = IdList()
feature_id_list.read(args.id_file)

with open(args.input, "r") as in_fd:
    with open(args.output, "w") as out_fd:
        for line in in_fd:
            if (line[0] == "#") and args.write_comments:
                out_fd.write(line)
                continue
            description_list = line.split("\t")[9].split(";")
            feature_id = description_list[0].split("=")[1]
            if feature_id not in feature_id_list:
                continue
            out_fd.write(line)
            while True:
                description_list = in_fd.next().split("\t")[9].split(";")




