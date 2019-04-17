#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Parsers.HGMD import HGMDvariants



parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input HGMD report",)
parser.add_argument("-v", "--variant_type", action="store", dest="variant_type", required=True,
                    help="Type of variant in hgmd report. Allowed: snp, indel, hgmd_indel")
parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator to use in output file. Default: tab")
args = parser.parse_args()

hgmd_variants = HGMDvariants(args.input, header=True)
coordinates_df = hgmd_variants.extract_variant_coordinates(args.variant_type)
coordinates_df.to_csv(args.output, sep=args.separator, index=False)

