#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse, os

from Bio import SeqIO
from Bio.SeqUtils import GC

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input file. Default - fasta")

args = parser.parse_args()

sequence_dict = SeqIO.index_db("temp.idx", args.input_file, args.format)

GC_content = [GC(sequence_dict[record].seq) for record in sequence_dict]


plt.figure(1, figsize=(6, 6))
plt.subplot(1, 1, 1)
plt.hist(GC_content)
plt.xlabel("GC content")
plt.ylabel("N")
plt.title("GC content")

for ext in ".png", ".svg", ".eps":
    plt.savefig(args.output_prefix + ext)

os.remove("temp.idx")