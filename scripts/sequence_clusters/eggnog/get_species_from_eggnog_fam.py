#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import shutil
import argparse

from Routines import EggNOGRoutines, FileRoutines
from CustomCollections.GeneralCollections import SynDict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Eggnog tsv file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Directory to write fam files named by species names")
parser.add_argument("-e", "--email", action="store", dest="email",
                    help="E-mail for request ot NCBI Taxonomy. If not set species ids will not be converted to species"
                         "names by queries to NCBI Taxonomy database(requires Internet connection)")

args = parser.parse_args()

EggNOGRoutines.get_species_from_eggnog_tsv(args.input, args.output_prefix, email=args.email)