#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Tools.WGA import LAST

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tab_file", action="store", dest="input_tab_file", required=True,
                    help="Input TAB file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output figure with dotplot")

parser.add_argument("-x", "--xsize", action="store", dest="xsize", default=2000, type=int,
                    help="X size of figure in pixels. Default: 2000")
parser.add_argument("-y", "--ysize", action="store", dest="ysize", default=2000, type=int,
                    help="Y size of figure in pixels. Default: 2000")

parser.add_argument("--seq_id_1", action="store", dest="seq_id_1",
                    type=LAST.split_string_by_comma,
                    help="Comma-separated list of sequence ids from genome one to show")
parser.add_argument("--seq_id_2", action="store", dest="seq_id_2",
                    type=LAST.split_string_by_comma,
                    help="Comma-separated list of sequence ids from genome two to show")
parser.add_argument("--seq_idfile_1", action="store", dest="seq_idfile_1",
                    help="File with sequence ids from genome one to show")
parser.add_argument("--seq_idfile_2", action="store", dest="seq_idfile_2",
                    help="File with sequence ids from genome two to show")
parser.add_argument("--seq_sort_1", action="store", dest="seq_sort_1", default="name",
                    help="Sequence sorting for genome one. Allowed: input, name(default), length, alignment")
parser.add_argument("--seq_sort_2", action="store", dest="seq_sort_2", default="name",
                    help="Sequence sorting for genome two. Allowed: input, name(default), length, alignment")


args = parser.parse_args()


LAST.plot(args.input_tab_file, args.output,
          first_genome_seq_id_list=args.seq_id_1 if args.seq_id_1 else IdList(filename=args.seq_idfile_1),
          second_genome_seq_id_list=args.seq_id_2 if args.seq_id_2 else IdList(filename=args.seq_idfile_2),
          first_genome_seq_order=args.seq_sort_1,
          second_genome_seq_order=args.seq_sort_2,
          xsize=args.xsize, ysize=args.ysize)
