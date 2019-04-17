#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Tools.Filter import Trimmomatic
from RouToolPa.Collections.General import TwoLvlDict
from RouToolPa.Routines.File import check_path




parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with sample reports")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples.")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write statistics")
parser.add_argument("-l", "--log_file", action="store", dest="log_file", default="trimmomatic.log",
                    help="Name of files with trimmomatic log. Default - trimmomatic.log")

args = parser.parse_args()

samples = sorted(args.samples.split(",") if args.samples else os.listdir(args.samples_dir))
present_samples = []
for sample in samples:
    if os.path.isdir(args.samples_dir + sample):
        present_samples.append(sample)

reports_dict = TwoLvlDict()

for sample in present_samples:
    print("Handling report from %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)
    trimmomatic_log = "%s/trimmomatic.log" % sample_dir
    reports_dict[sample] = Trimmomatic.parse_log(trimmomatic_log)

reports_dict.write(args.output)
