#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import IdList



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--fam_file", action="store", dest="fam_file", required=True,
                    help="File with families")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write ids")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

id_list = IdList()
id_list.read(args.fam_file, close_after_if_file_object=True, column_number=1, id_in_column_separator=",")
id_list.write(args.output, close_after_if_file_object=True)
