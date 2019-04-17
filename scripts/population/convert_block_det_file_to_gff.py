#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Population import PLINK

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--block_det_file", action="store", dest="block_det_file", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of PLINK .block.det files")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output GFF file with block")
parser.add_argument("-f", "--output_format", action="store", dest="output_format", default="gff",
                    help="Format of output file. Default: gff, bed")
parser.add_argument("-b", "--haploblock_id_prefix", action="store", dest="haploblock_id_prefix", default="HB",
                    help="Prefix of haploblock ids. Default: HB")

args = parser.parse_args()

PLINK.convert_blocks_det_file_to_gff(args.block_det_file, args.output, output_format=args.output_format,
                                     haploblock_id_prefix=args.haploblock_id_prefix)
