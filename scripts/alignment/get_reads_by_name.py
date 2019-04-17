#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Tools.Samtools import SamtoolsV1
from RouToolPa.Collections.General import IdList



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Input sam file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with reads. Default: stdout")
parser.add_argument("-r", "--read_name_file", action="store", dest="read_name_file", required=True,
                    help="File with full read names or their fragments")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="include",
                    help="Output mode. Allowed: include(default), remove")
parser.add_argument("-c", "--comparison_mode", action="store", dest="comparison_mode",
                    default="exact",
                    help="Read name comparison mode. Allowed: exact(default), partial")

args = parser.parse_args()

input_sam_fd = open(args.input, "r") if args.input else sys.stdin
output_sam_fd = open(args.output, "w") if args.output else sys.stdout

read_name_list = IdList(filename=args.read_name_file)
SamtoolsV1.get_reads_by_name(read_name_list, input_sam_fd, output_sam_fd,
                             mode=args.mode, search_mode=args.comparison_mode)

if args.input:
    input_sam_fd.close()
if args.output:
    output_sam_fd.close()
