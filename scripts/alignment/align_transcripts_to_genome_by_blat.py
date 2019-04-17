#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Alignment import BLAT

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--database", action="store", dest="database", required=True,
                    help="Input file with database")
parser.add_argument("-q", "--query_fasta", action="store", dest="query_fasta", required=True,
                    help="Fasta file with query sequences")

parser.add_argument("-p", "--path_to_blat_dir", action="store", dest="path_to_blat_dir",
                    default="", help="Path to BLAT directory")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=4,
                    help="Number of threads")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

BLAT.threads = args.threads
BLAT.path = args.path_to_blat_dir
BLAT.align_transcripts_for_scaffolding(args.database, args.query_fasta, args.output, split_dir="splited_input/",
                                       splited_output_dir="splited_output_dir/", threads=None, remove_tmp_dirs=True,
                                       async_run=False, external_process_pool=None)
