#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import FastQRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_list", action="store", dest="forward_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with forward sequences.")
parser.add_argument("-r", "--reverse_list", action="store", dest="reverse_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with forward sequences.")
parser.add_argument("-i", "--index_list", action="store", dest="index_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated ist of files with index sequences.")
parser.add_argument("-b", "--barcode_file", action="store", dest="barcode_file", required=True,
                    help="File with 10X barcodes")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

FastQRoutines.extract_10x_barcodes(args.forward_list, args.reverse_list, args.index_list, args.barcode_file,
                                   args.output_prefix, buffering=100000000, read_index_length=16, linker_length=6,
                                   min_forward_read_len=50)
