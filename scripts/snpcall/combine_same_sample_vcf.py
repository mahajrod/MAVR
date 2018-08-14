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


args = parser.parse_args()

VCFRoutines.combine_same_samples_vcfs(args.vcf_list, args.output, close_fd_after=False, extension_list=[".vcf", ])
