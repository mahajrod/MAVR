#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.Population import PLINK

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--block_det_file", action="store", dest="block_det_file", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of PLINK .block.det files")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output GFF file with block")

args = parser.parse_args()

PLINK.convert_blocks_det_file_to_gff(args.block_det_file, args.output)
