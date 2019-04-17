#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines.File import make_list_of_path_to_files



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
parser.add_argument("-c", "--case_insensitive", action="store_true", dest="case_insensitive",
                    help="Case insensitive comparison. With this flag all ids are converted to upper case "
                         "before comparison and are reported in upper case too. Default - False")

args = parser.parse_args()


Tool.intersect_ids_from_files(args.group_a_list, args.group_b_list, args.output, mode=args.mode,
                              case_insensitive=args.case_insensitive)
