#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Routines import VCFRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")
parser.add_argument("-i", "--vcf_list", action="store", dest="vcf_list", required=True,
                    type=VCFRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of vcf files")
parser.add_argument("-s", "--sort", action="store_true", dest="sort", default=False,
                    help="Sort vcf files. Default:False")

args = parser.parse_args()

VCFRoutines.combine_same_samples_vcfs(args.output, vcf_list=args.vcf_list, close_fd_after=False, extension_list=[".vcf", ],
                                      order_vcf_files=False, sort=args.sort)
