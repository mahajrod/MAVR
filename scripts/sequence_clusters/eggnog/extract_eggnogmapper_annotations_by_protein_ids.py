#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Emapper


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with emapper annotations")
parser.add_argument("-d", "--protein_id_file", action="store", dest="protein_id_file", required=True,
                    help="File with protein ids to extract")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

Emapper.extract_eggnogmapper_annotations_by_protein_ids(args.input, args.protein_id_file, args.output)
