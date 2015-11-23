#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines.File import make_list_of_path_to_files
from CustomCollections.GeneralCollections import IdSet


parser = argparse.ArgumentParser()

parser.add_argument("-w", "--white_list", action="store", dest="white_list", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with ids from white list")
parser.add_argument("-b", "--black_list", action="store", dest="black_list", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with ids from black list")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with allowed ids")

args = parser.parse_args()

white_lists = []
black_lists = []

white_set = set()
black_set = set()

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

