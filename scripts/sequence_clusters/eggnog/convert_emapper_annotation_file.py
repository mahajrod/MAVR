#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Emapper


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with emapper annotations")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", required=True,
                    help="Prefix of output files")
parser.add_argument("-p", "--eggnogdb_prefix", action="store", dest="eggnogdb_prefix",
                    help="EggNOG database prefix for clusters. Default: None")
parser.add_argument("-s", "--species_name", action="store", dest="species_name",
                    help="Species name to use for labeling of proteins. Default: not set")


args = parser.parse_args()

Emapper.convert_emapper_annotation_file(args.input, args.output, eggnogdb_prefix=args.eggnogdb_prefix,
                                        species_name=args.species_name, label_separator=".")
