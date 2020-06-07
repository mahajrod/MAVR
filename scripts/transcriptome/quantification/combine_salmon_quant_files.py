#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Transcriptome import Salmon


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file_list", action="store", dest="file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with Salmon quants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output_files")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample names."
                         "If not set filenames will be treated as sample names")
parser.add_argument("-t", "--tpm_threshold_list", action="store", dest="tpm_threshold_list", default=[1.0, 0.5, 0.1],
                    type=lambda s: list(map(float, s.split(","))),
                    help="List of thresholds for TPM. Default: 1.0, 0.5, 0.1")

args = parser.parse_args()

Salmon.combine_per_sample_counts(args.file_list,
                                 args.sample_list if args.sample_list else args.file_list,
                                 args.output_prefix, tpm_threshold_list=args.tpm_threshold_list)
