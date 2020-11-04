#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import HaplotypeRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--alignment", action="store", dest="aln_file", required=True,
                    help="Input file with alignment in fasta format")
parser.add_argument("-y", "--hap_file", action="store", dest="hap_file",
                    help="Fam file with haplotypes")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-t", "--traits_file", action="store", dest="traits_file",
                    help="Tab-separated file with traits for each haplotype")
parser.add_argument("-w", "--whitelist_file", action="store", dest="whitelist_file",
                    help="File with ids from whitelist. Default: not set")

args = parser.parse_args()

HaplotypeRoutines.prepare_template_for_popart(args.aln_file,
                                              args.output,
                                              haplotype_fam_file=args.hap_file,
                                              traits_file=args.traits_file,
                                              whitelist_file=args.whitelist_file)
