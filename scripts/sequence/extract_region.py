#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from RouToolPa.Routines.File import make_list_of_path_to_files



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
        regions_list[i].append(None)
    return regions_list


def get_regions_from_file(filename):
    regions_list = []

    with open(filename, "r") as region_fd:
        for line in region_fd:
            tmp = line.strip().split("\t")
            regions_list.append((tmp[0], (int(tmp[1]), int(tmp[2])),
                                 True if tmp[3] == "-" else None,
                                 tmp[4] if len(tmp) >= 5 else None))

    return regions_list


def region_generator(sequence_dict, region_list):
    for record_id, (start, stop), reverse, region_name in region_list:
        if record_id not in sequence_dict:
            if args.output != "stdout":
                print("Record %s is absent" % record_id)
            continue
        seq = sequence_dict[record_id].seq[start-1: stop]
        seq_id = "%s_%i-%i" % (record_id, start, stop)
        if reverse:
            seq = seq.reverse_complement()
            seq_id += "_r"
        if region_name:
            seq_id += "_%s" % region_name

        yield SeqRecord(seq=seq, id=seq_id)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with sequence of region")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed formats genbank, fasta(default)")
parser.add_argument("-r", "--regions", action="store", dest="regions",
                    type=get_regions_from_string,
                    help="Comma-separated_list of regions to extract. Option is ignored if -l/--region_file "
                         "option is used"
                         "Format of region: <record_id>:<start>-<stop>[:r]. Coordinates MUST be one-based. "
                         "Place 'r' after coordinates if region should be reverse complement ")
parser.add_argument("-l", "--region_file", action="store", dest="region_file",
                    type=get_regions_from_file,
                    help="Bed-like file with regions to extract. This option overrides -r/--regions option")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

if (not args.regions) and (not args.region_file):
    raise ValueError("Neither -r/--regions nor -l/--region_file options were set")

args.input = make_list_of_path_to_files(args.input)
tmp_index_file = "temp.idx"

sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)

SeqIO.write(region_generator(sequence_dict, args.regions if args.regions else args.region_file),
            out_fd, format=args.format)
os.remove(tmp_index_file)

if args.output != "stdout":
    out_fd.close()

