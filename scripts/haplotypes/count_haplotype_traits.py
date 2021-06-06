#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import HaplotypeRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-y", "--hap_file", action="store", dest="hap_file", required=True,
                    help="File with haplotype traits")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")


args = parser.parse_args()

HaplotypeRoutines.count_haplotype_traits(args.hap_file, args.output)
