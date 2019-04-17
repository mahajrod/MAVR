#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from BCBio import GFF



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with transcript lengths")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

for record in GFF.parse(open(args.input_gff)):
    for feature in record.features:
        #print feature
        if feature.type == "gene":
            for subfeature in feature.sub_features:
                #print subfeature
                print(subfeature.type)
                if subfeature.type == "mRNA" or subfeature.type == "transcript":
                    out_fd.write("%s\t%s\t%i\n" % (feature.id, subfeature.id, len(subfeature)))

out_fd.close()



