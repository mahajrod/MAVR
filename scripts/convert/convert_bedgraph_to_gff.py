#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input bedgraph file(0-based, python notation)")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gff file")
parser.add_argument("-f", "--feature_type", action="store", dest="feature_type", required=True,
                    help="Feature type to use in gff file")
parser.add_argument("-s", "--source", action="store", dest="source", default="source",
                    help="Source to use in gff file")
parser.add_argument("-d", "--id_prefix", action="store", dest="id_prefix",
                    default="ID", help="Id prefix for gff file")

args = parser.parse_args()

AnnotationsRoutines.convert_bedgraph_to_gff(args.input, args.output, args.feature_type,
                                            id_prefix=args.id_prefix, source=args.source)
