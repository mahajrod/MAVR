#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chunks_dir", action="store", dest="chunks_dir", required=True,
                    help="Directory with chunks")
parser.add_argument("-p", "--chunks_prefix", action="store", dest="chunks_prefix",
                    help="Prefix of chunk files")
parser.add_argument("-u", "--chunks_suffix", action="store", dest="chunks_suffix",
                    help="Suffix of chunk files")

parser.add_argument("-a", "--starting_chunks_number", action="store", dest="starting_chunks_number", type=int,
                    help="Starting chunk number")
parser.add_argument("-b", "--ending_chunks_number", action="store", dest="ending_chunks_number", type=int,
                    help="Ending chunk number")
parser.add_argument("-n", "--chunks_number_list", action="store", dest="chunks_number_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of  chunk numbers")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output merged file")
parser.add_argument("-e", "--header_prefix", action="store", dest="header_prefix", default="#",
                    help="Header prefix")
parser.add_argument("-s", "--sorting_options", action="store", dest="sorting_options",
                    help="Sorting options for sort utility. Default: not set, i.e. no sort")
parser.add_argument("-r", "--separator", action="store", dest="separator", default="_",
                    help="Separator between prefix and chunk number in chunk filename. Default: '_'")

args = parser.parse_args()

FileRoutines.combine_chunks_with_header(args.chunks_dir,
                                        args.chunks_prefix,
                                        args.output,
                                        starting_chunk=args.starting_chunks_number,
                                        end_chunk=args.ending_chunks_number,
                                        chunk_number_list=args.chunks_number_list,
                                        chunks_suffix=args.chunks_suffix,
                                        header_prefix=args.header_prefix,
                                        sorting_options=args.sorting_options,
                                        separator=args.separator)
