#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.RepeatMasking import RepeatModeler

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file/directory with files in fasta format. Only files with extensions"
                         " .fa/.fasta are recognized.")
parser.add_argument("-n", "--database_name", action="store", dest="database_name", required=True,
                    help="Name of database to use")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")

parser.add_argument("-e", "--engine", action="store", dest="engine", default="ncbi",
                    help="Engine to use for repeat annotation. Default - 'ncbi'")
parser.add_argument("-f", "--repeat_families", action="store", dest="repeat_families", required=True,
                    help="File to write annotated repeat families")
args = parser.parse_args()

RepeatModeler.threads = args.threads
RepeatModeler.build_db(args.database_name, fasta_dir=args.input, file_with_filenames=None, engine=args.engine)
RepeatModeler.annotate_repeats(args.database_name, engine=args.engine, recover_dir=None)
