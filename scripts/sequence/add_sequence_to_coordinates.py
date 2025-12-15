#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os
import sys

import pandas as pd
from Bio import SeqIO
from BCBio import GFF


from RouToolPa.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input_fasta", default=sys.stdin,
                    help="Input fasta file with sequences. Default: STDIN")
parser.add_argument("-c", "--coordinate_file", action="store", dest="coordinate_file", required=True,
                    help="Input file with coordinates")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Column separator in the file with coordinates. Default: TAB")
parser.add_argument("-b", "--one-based", action="store_true", dest="one_based",
                    help="Flag to indicate that coordinates are 1-based and regions are closed (end coordinate is a part of region). "
                         "By default coordinates are treated as ZERO-BASED and HALF-OPEN like in .bed format."
                         "Set this flag if, for example, your coordinates are from .vcf or .gff file")
parser.add_argument("--header",  action="store_true", dest="header",
                    help="Flag to indicate that file with coordinates has a header")
parser.add_argument("--scaffold_id_column_idx",  action="store", dest="scaffold_id_column_idx", default=0, type=int,
                    help="Index (zero-based) of column containing scaffold ids in coordinate file. Default: 0, like in .bed files")
parser.add_argument("--start_column_idx",  action="store", dest="start_column_idx", default=1, type=int,
                    help="Index (zero-based) of column containing start coordinate of regions in coordinate file. "
                         "Default: 1, like in .bed files")
parser.add_argument("--end_column_idx",  action="store", dest="end_column_idx", default=2, type=int,
                    help="Index (zero-based) of column containing end coordinate of regions in coordinate file. "
                         "Default: 2, like in .bed files")
parser.add_argument("--strand_column_idx",  action="store", dest="strand_column_idx", default=None, type=int,
                    help="Index (zero-based) of column containing strand of regions in coordinate file. "
                         "Default: not set. If not set, sequences will be extracted from '+' strand, "
                         "but 'strand' value in the output file will be set to NA")
parser.add_argument("--region_id_column_idx",  action="store", dest="region_id_column_idx", default=None, type=int,
                    help="Index (zero-based) of column containing ids of regions in coordinate file. "
                         "Default: not set. If not set, regions will get consequent ids starting from 1 with prefix 'LOC', "
                         "i.e. LOC1, LOC2, ..., LOC100500, etc")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: STDOUT")


args = parser.parse_args()



fasta_collection = CollectionSequence(in_file=args.input_fasta, format="fasta")

columns_to_read = [args.scaffold_id_column_idx, args.start_column_idx, args.end_column_idx]
if args.strand_column_idx is not None:
    columns_to_read.append(args.strand_column_idx)
if args.region_id_column_idx is not None:
    columns_to_read.append(args.region_id_column_idx)
    region_id_column_indexxxxxx = 3 if args.strand_column_idx is None else 4
else:
    region_id_column_indexxxxxx = None

coordinate_df = pd.read_csv(args.coordinate_file, sep=args.separator, header=0 if args.header else None,
                            usecols=columns_to_read)

if args.one_based:
    coordinate_df.iloc[:, 1] = coordinate_df.iloc[:, 1] - 1

output_df = fasta_collection.get_sequence_for_regions(coordinate_df,
                                                     scaffold_id_column_idx=0,
                                                     start_column_idx=1,
                                                     end_column_idx=2,
                                                     strand_column_idx=3 if args.strand_column_idx is not None else None,
                                                     region_id_column_idx=region_id_column_indexxxxxx)

output_df.to_csv(args.output, sep="\t", header=True, index=False, na_rep=".")