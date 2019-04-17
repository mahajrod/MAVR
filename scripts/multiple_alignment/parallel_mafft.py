#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.MultipleAlignment import MAFFT
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", type=lambda x: x.split(","),
                    help="Comma-separated list of files or directory with files containing sequences to be aligned")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads per alignment")
parser.add_argument("-p", "--processes", action="store", dest="processes", type=int, default=1,
                    help="Number of simalteneously running alignments")
parser.add_argument("-o", "--output_directory", action="store", dest="output", type=FileRoutines.check_path,
                    help="Output directory")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="globalpair",
                    help="Alignment mode. Default: 'globalpair'. Allowed: globalpair, localpair, genafpair")
parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                    help="Quiet output")
parser.add_argument("-x", "--maxiterate", action="store", dest="maxiterate", type=float,
                    help="Maximum number of iterations")
parser.add_argument("-f", "--offset", action="store", dest="offset", type=float,
                    help="Offset (works like gap extension penalty)")
parser.add_argument("-g", "--gap_open_penalty", action="store", dest="gap_open_penalty", type=float,
                    help="Gap open penalty")

args = parser.parse_args()

FileRoutines.safe_mkdir(args.output)

MAFFT.threads = args.threads
MAFFT.parallel_align(FileRoutines.make_list_of_path_to_files(args.input), args.output, output_suffix="alignment",
                     gap_open_penalty=args.gap_open_penalty, offset=args.offset, maxiterate=args.maxiterate,
                     quiet=args.quiet, mode=args.mode, number_of_processes=args.processes, anysymbol=True)
