#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from BConverters import SequenceConverters

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input sequence file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output sequence file")
parser.add_argument("-f", "--input_format", action="store", dest="input_format", required=True,
                    help="Format of input sequence file")
parser.add_argument("-g", "--output_format", action="store", dest="output_format", required=True,
                    help="Format of output sequence file")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode",
                    default="parse", help="Parsing mode. Allowed: parse(default), index, index_db")
parser.add_argument("-a", "--alphabet", action="store", dest="alphabet",
                    help="Sequence alphabet. Default: not set")
args = parser.parse_args()

SequenceConverters.convert_sequences(args.input, args.input_format, args.output, args.output_format,
                                     parsing_mode=args.parsing_mode, input_index="tmp.idx", alphabet=args.alphabet)
