#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines import SequenceRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta file with reference")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-s", "--split_scaffolds", action="store_true", dest="split_scaffolds", default=False,
                    help="Split scaffolds. Default: False")
parser.add_argument("-m", "--max_length", action="store", dest="max_length", type=int,
                    help="Soft maximum length of region(1.5x longer regions are allowed). Default: not set")
parser.add_argument("-n", "--max_seq_number", action="store", dest="max_seq_number", type=int, default=1,
                    help="Maximum number of sequences per region. Default: 1")
parser.add_argument("-b", "--scaffold_black_list_file", action="store", dest="scaffold_black_list_file",
                    type=lambda s: IdList(filename=s),
                    help="File with scaffolds from black list")
parser.add_argument("-x", "--min_scaffold_len", action="store", dest="min_scaffold_len", type=int, default=None,
                    help="Minimum length of scaffold to be included in regions. Default: not set")

args = parser.parse_args()


SequenceRoutines.prepare_region_list_by_length(max_length=args.max_length,
                                               max_seq_number=args.max_seq_number,
                                               length_dict=None,
                                               reference=args.reference,
                                               parsing_mode="parse",
                                               output_dir=args.output_dir,
                                               split_scaffolds=args.split_scaffolds,
                                               min_scaffold_length=args.min_scaffold_len,
                                               black_list_scaffolds=None)


