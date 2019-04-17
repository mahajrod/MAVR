#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from Bio import SeqIO
from RouToolPa.Tools.Bedtools import Flank, GetFasta
from RouToolPa.Routines.Sequence import SequenceRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fasta", action="store", dest="input_fasta",
                    help="fasta file with sequences")
parser.add_argument("-b", "--bed", action="store", dest="bed",
                    help="bed file with features")
parser.add_argument("-o", "--output", action="store", dest="output_fasta",
                    default="out.fasta", help="output file - default: out.fasta.")
parser.add_argument("-l", "--left", action="store", dest="left",
                    help="The number of base pairs that a flank should start from orig. "
                         "start coordinate.- (Integer) or (Float, e.g. 0.1) if used with -pct.")
parser.add_argument("-r", "--right", action="store", dest="right",
                    help="The number of base pairs that a flank should end from orig. end coordinate")
parser.add_argument("--pct", action="store_true", dest="fraction_mode",
                    help="Define -l and -r as a fraction of the feature's length."
                         "E.g. if used on a 1000bp feature, -l 0.50,will add 500 bp \"upstream\".  Default = false")
parser.add_argument("-s", action="store_true", dest="use_strand",
                    help="Define -l and -r based on strand. E.g. if used, -l 500 for a "
                         "negative-stranded feature, it will start the flank 500 bp downstream.  Default = false.")

args = parser.parse_args()

record_dict = SeqIO.index_db("temp_index.idx", [args.input_fasta], format="fasta")

SequenceRoutines.get_lengths(record_dict, out_file="fasta_lengths.t", write=True)

if args.fraction_mode:
    left = float(args.left)
    right = float(args.right)
else:
    left = int(args.left)
    right = int(args.right)

Flank.get(args.bed, "fasta_lengths.t", left, right, fraction_mode=args.fraction_mode,
          strand_based=args.use_strand, out_file="flank.t")

GetFasta.get("flank.t", args.input_fasta, args.output_fasta, use_strand=args.use_strand)
