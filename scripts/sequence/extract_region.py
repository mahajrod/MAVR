#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Routines.File import read_ids, make_list_of_path_to_files
from Routines.Sequence import record_by_id_generator


def get_regions_from_string(s):
    regions_list = s.split(",")
    for i in range(0, len(regions_list)):
        regions_list[i] = regions_list[i].split(":")
        if len(regions_list[i]) > 2:
            if regions_list[i][2] != "r":
                raise ValueError("Error in region")
        elif len(regions_list[i]) == 2:
                regions_list[i].append(None)
        regions_list[i][1] = map(int, regions_list[i][1].split("-"))
    return regions_list


def region_generator(sequence_dict, region_list):
    for record_id, (start, stop), reverse in region_list:
        if record_id not in sequence_dict:
            if args.output != "stdout":
                print("Record %s is absent" % record_id)
            continue
        seq = sequence_dict[record_id].seq[start-1: stop]
        seq_id = "%s_%i-%i" % (record_id, start, stop)
        if reverse:
            seq = seq.reverse_complement()
            seq_id += "_r"

        yield SeqRecord(seq=seq, id=seq_id)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with sequence of region")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed formats genbank, fasta(default)")
parser.add_argument("-r", "--regions", action="store", dest="regions", required=True,
                    type=get_regions_from_string,
                    help="Comma-separated_list of regions to extract. "
                         "Format of region: <record_id>:<start>-<stop>[:r]. Coordinates MUST be one-based. "
                         "Place 'r' after coordinates if region should be reverse complement ")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

args.input = make_list_of_path_to_files(args.input)
tmp_index_file = "temp.idx"

sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)

SeqIO.write(region_generator(sequence_dict, args.regions), out_fd, format=args.format)
os.remove(tmp_index_file)

if args.output != "stdout":
    out_fd.close()

