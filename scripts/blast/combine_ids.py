#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.Abstract import Tool
from Routines.File import make_list_of_path_to_files
from CustomCollections.GeneralCollections import IdSet


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--group_a_list", action="store", dest="group_a_list", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with ids from group A")
parser.add_argument("-b", "--group_b_list", action="store", dest="group_b_list", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with ids from group B")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with allowed ids. Default - stdout")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="common",
                    help="Operation mode. Possible variants: common, only_a, only_b, not_common, combine, count"
                         ". Default - common")

args = parser.parse_args()

#white_lists = []
#black_lists = []

Tool.intersect_ids_from_files(args.white_list, args.black_list, args.output, mode=args.mode)

"""
white_set = IdSet()
black_set = IdSet()

for filename in args.white_list:
    id_set = IdSet()
    id_set.read(filename, comments_prefix="#")
    white_set = white_set | id_set

for filename in args.black_list:
    id_set = IdSet()
    id_set.read(filename, comments_prefix="#")
    black_set = black_set | id_set

final_set = IdSet(white_set - black_set)

final_set.write(args.output)
"""
