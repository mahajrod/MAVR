#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Exonerate
from RouToolPa.Routines.File import make_list_of_path_to_files


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Input comma-separated list of files/directories with exonerate output")
parser.add_argument("-r", "--reference_protein_file", action="store", dest="reference_protein_file", required=True,
                    help="File with reference proteins used for alignment")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-g", "--gene_id_prefix", action="store", dest="gene_id_prefix", default="GEN",
                    help="Prefix of gene id to use. Default: GEN")
parser.add_argument("-t", "--transcript_id_prefix", action="store", dest="transcript_id_prefix", default="TR",
                    help="Prefix of transcript id to use. Default: TR")
parser.add_argument("-n", "--number_digits_in_id", action="store", dest="number_digits_in_id", default=8, type=int,
                    help="Number of digits in id. Default: 8")
parser.add_argument("-q", "--query_gff_presence", action="store_true", dest="query_gff_presence",
                    help="Presence of query gff in input file. Default: False")
args = parser.parse_args()

Exonerate.split_output(args.input, args.output_prefix, args.reference_protein_file, gene_prefix=args.gene_id_prefix,
                       transcript_prefix=args.gene_id_prefix, number_len=args.number_digits_in_id,
                       query_gff_presence=args.query_gff_presence)
