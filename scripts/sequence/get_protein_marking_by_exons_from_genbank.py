#!/usr/bin/env python
__author__ = 'mahajrod'
import argparse
from Bio import SeqIO
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Genbank file with annotations")
parser.add_argument("--fast_parsing", action="store_true", dest="fast_parsing",
                    help="Fast parsing mode - high memory consumption. Default: false")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    default="output", help="Prefix of output files")

args = parser.parse_args()


record_dict = SeqIO.to_dict(SeqIO.parse(args.input, format="genbank")) if args.fast_parsing else SeqIO.index_db("temp_index.idx", [args.input], format="genbank")

SequenceRoutines.get_protein_marking_by_exons_from_genbank(record_dict, args.output_prefix,
                                                           protein_id_field_in_cds_feature="protein_id")

#os.remove("temp_index.idx")
