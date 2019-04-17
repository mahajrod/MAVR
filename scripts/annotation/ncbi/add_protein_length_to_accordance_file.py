#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import SynDict


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--accordance_file", action="store", dest="accordance_file", required=True,
                    help="Gene to protein accordance file")
parser.add_argument("-l", "--len_file", action="store", dest="length_file", required=True,
                    help="File with protein lengths")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

length_dict = SynDict(filename=args.length_file)

with open(args.accordance_file, 'r') as in_fd:
    with open(args.output, 'w') as out_fd:
        for line in in_fd:
            line_list = line.strip().split('\t')
            protein_id = line_list[1]
            if protein_id not in length_dict:
                print("No length for protein %s" % protein_id)
                continue
            out_fd.write("%s\t%s\t%s\n" % (line_list[0], line_list[1], length_dict[protein_id]))
