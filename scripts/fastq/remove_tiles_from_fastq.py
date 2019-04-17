#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_reads", action="store", dest="forward_reads", required=True,
                    help="Forward reads")
parser.add_argument("-r", "--reverse_reads", action="store", dest="reverse_reads", required=True,
                    help="Reverse reads")
parser.add_argument("-a", "--forward_tiles_to_remove", action="store", dest="forward_tiles_to_remove",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of tiles to remove from forward reads")
parser.add_argument("-b", "--reverse_tiles_to_remove", action="store", dest="reverse_tiles_to_remove",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of tiles to remove from reverse reads")
parser.add_argument("-o", "--output-prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

if args.forward_tiles_to_remove is None:
    args.forward_tiles_to_remove = []

if args.reverse_tiles_to_remove is None:
    args.reverse_tiles_to_remove = []

if (not args.reverse_tiles_to_remove) and (not args.forward_tiles_to_remove):
    raise ValueError("Tile black lists were not set for both forward and reverse reads")

FastQRoutines.remove_tiles_from_fastq(args.forward_reads, args.forward_tiles_to_remove,
                                      args.reverse_reads, args.reverse_tiles_to_remove, args.output_prefix)
