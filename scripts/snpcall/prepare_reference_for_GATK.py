#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="File with reference genome")
parser.add_argument("-s", "--samtools_directory", action="store", dest="samtools_dir", default="",
                    help="Directory with samtools binaries")
parser.add_argument("-p", "--picard_directory", action="store", dest="picard_dir", default="",
                    help="Directory with PICARD jar")

args = parser.parse_args()


SequenceRoutines.prepare_reference_for_GATK(args.reference,
                                            picard_dir=args.picard_dir,
                                            samtools_dir=args.samtools_dir)
