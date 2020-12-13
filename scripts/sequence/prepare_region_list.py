#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta file with reference")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    help="Output directory. If not set output will be written to stdout")
parser.add_argument("-s", "--split_scaffolds", action="store_true", dest="split_scaffolds", default=False,
                    help="Split scaffolds. Default: False")
parser.add_argument("-m", "--max_length", action="store", dest="max_length", type=int,
                    help="Soft maximum length of region(1.5x longer regions are allowed). Default: not set")
parser.add_argument("-n", "--max_seq_number", action="store", dest="max_seq_number", type=int, default=1,
                    help="Maximum number of sequences per region. Default: 1")
parser.add_argument("-b", "--scaffold_black_list_file", action="store", dest="scaffold_black_list_file",
                    type=lambda s: IdList(filename=s) if os.path.isfile(s) else IdList(s.split(",")),
                    help="File or comma-separated list with scaffolds from black list")
parser.add_argument("-w", "--scaffold_white_list_file", action="store", dest="scaffold_white_list_file",
                    type=lambda s: IdList(filename=s) if os.path.isfile(s) else IdList(s.split(",")),
                    help="File or comma-separated list with scaffolds from white list")
parser.add_argument("-x", "--min_scaffold_len", action="store", dest="min_scaffold_len", type=int, default=None,
                    help="Minimum length of scaffold to be included in regions. Default: not set")
parser.add_argument("-g", "--region_file_format", action="store", dest="region_file_format", default='simple',
                    help="Output region file format. "
                         "Allowed: 'simple' (default, not appliable for stdout), 'GATK', 'samtools'")
args = parser.parse_args()


SequenceRoutines.prepare_region_list_by_length(max_length=args.max_length,
                                               max_seq_number=args.max_seq_number,
                                               length_dict=None,
                                               reference=args.reference,
                                               parsing_mode="parse",
                                               output_dir=args.output_dir,
                                               split_scaffolds=args.split_scaffolds,
                                               min_scaffold_length=args.min_scaffold_len,
                                               black_list_scaffolds=args.scaffold_black_list_file,
                                               white_list_scaffolds=args.scaffold_white_list_file,
                                               region_file_format=args.region_file_format)


