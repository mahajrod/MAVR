#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BCBio import GFF


from Routines.Sequence import record_generator


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with extracted_annotations")

parser.add_argument("-t", "--annotation_types", action="store", dest="types", default="gene",
                    type=lambda s: s.split(),
                    help="Comma-separated list of annotation types to extract")

args = parser.parse_args()

for record in GFF.parse(open(args.input_gff)):
    print record.id



