#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from CustomCollections.GeneralCollections import TwoLvlDict
from collections import OrderedDict

from Tools.Filter import Trimmomatic

from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with sample reports")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples.")
parser.add_argument("-o", "--output", action="store", dest="output_dir", required=True,
                    help="File to write statistics")
parser.add_argument("-l", "--log_file", action="store", dest="log_file", default="trimmomatic.log",
                    help="Name of files with trimmomatic log")
args = parser.parse_args()

Trimmomatic.jar_path = args.path_to_trimmomatic_dir

samples = args.samples.split(",") if args.samples else os.listdir(args.samples_dir)
reports_dict = TwoLvlDict()

for sample in samples:
    print("Handling report from %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)
    trimmomatic_log = "%s/trimmomatic.log" % sample_dir
    reports_dict[sample] = Trimmomatic.parse_log(trimmomatic_log)

reports_dict.write(args.output_dir)
