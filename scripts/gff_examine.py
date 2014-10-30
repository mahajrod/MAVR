#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
import pprint
from BCBio.GFF import GFFExaminer

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--gff", action="store", dest="gff",
                    help="gff to examine")
args = parser.parse_args()

examiner = GFFExaminer()

with open(args.gff, "r") as in_fd:
    pprint.pprint(examiner.parent_child_map(in_fd))