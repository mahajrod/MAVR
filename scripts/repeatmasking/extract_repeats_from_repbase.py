#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species", action="store", dest="species",
                    help="The full name ( case insensitive ) of the species you would like to"
                          "search for in the database. This will return all the repeats which"
                          "would be used in a RepeatMasker search against this species. This"
                          "includes repeats contained in the clade given by 'species name' and"
                          "ancestral repeats of 'species name'. Lastly ubiquitous sequences"
                          "such as RNAs and simple repeats are also included.")

parser.add_argument("-c", "--clade", action="store", dest="clade",
                    help="This will modify the default behaviour of the species option and"
                         "return only the repeats which are specific to your species or any of"
                         "it descendents. This is useful for identifying how rich the database"
                         "of repeats is for a given species/clade.")

parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

RepeatMasker.extract_repeats_from_database(args.output, species=args.species, clade=args.clade, stat_mode=None)
