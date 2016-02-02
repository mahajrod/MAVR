#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from CustomCollections.GeneralCollections import SynDict, IdSet


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fam file")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids of families to extract")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write extracted families. Default - stdout")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print not found ids. Default - no")

args = parser.parse_args()

out_file = sys.stdout if args.output == "stdout" else open(args.output, "w")

fam_dict = SynDict()
fam_dict.read(args.input)

id_set = IdSet()
id_set.read(args.id_file)

extracted_dict = SynDict()
for id_entry in id_set:
    if id_entry in fam_dict:
        extracted_dict[id_entry] = fam_dict[id_entry]
    else:
        if args.verbose:
            print("%s was not found" % id_entry)

extracted_dict.write(out_file, close_after_if_file_object=True)




