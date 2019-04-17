#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Tools.HMMER import HMMER3
from RouToolPa.Collections.General import IdList



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with hits")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids")
parser.add_argument("-e", "--header", action="store_true", dest="header",
                    help="Header is present in id file. Default: False")
parser.add_argument("-f", "--format", action="store", dest="format", required=True,
                    help="Format of the file with hits")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
id_list = IdList()
id_list = id_list.read(args.id_file, header=args.header)

HMMER3.extract_hits_by_query_ids(id_list, args.input, args.output,
                                 fileformat=args.format,
                                 close_after_if_file_object=True)

out_fd.close()
