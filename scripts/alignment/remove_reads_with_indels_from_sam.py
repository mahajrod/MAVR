#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
import sys
import pysam


parser = argparse.ArgumentParser()

parser.add_argument("-s", "--input_sam", action="store", dest="input_sam", default=sys.stdin,
                    help="Input sam file. Default: stdin.")
parser.add_argument("-o", "--output_sam", action="store", dest="output_sam", default=sys.stdout,
                    help="Output sam file. Default: stdout.")
parser.add_argument("-x", "-max_length", action="store", dest="max_length", default=100,  type=int,
                    help="Maximal length of indel. Default: 100")
#parser.add_argument("-m", "-min_mapq", action="store", dest="min_mapq", default=0,  type=int,
#                    help="Minimum mappinq quality required to retain the alignment. "
#                         "Default: 0, i.e all reads will be retained")


args = parser.parse_args()

in_samfile = pysam.AlignmentFile(args.input_sam, "r")
out_samfile = pysam.AlignmentFile(args.output_sam, "w", template=in_samfile)

for aln in in_samfile.fetch():
    if aln.cigartuples is None:
        continue
    for operation, length in aln.cigartuples:
        if ((operation == 1) or (operation == 2)) and (length > args.max_length):
            break
    else:
        out_samfile.write(aln)

out_samfile.close()
in_samfile.close()
