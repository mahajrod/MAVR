#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import VCFRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output vcf. Default: stdout")

parser.add_argument("-d", "--id_prefix", action="store", dest="id_prefix", required=True,
                    help="Prefix for variant ids")
parser.add_argument("-r", "--retain_old_id", action="store_true", dest="retain_old_id", default=True,
                    help="Retain old ids of variants. Default: True")

args = parser.parse_args()

VCFRoutines.add_variant_ids_to_vcf(args.input, args.output, id_prefix=args.id_prefix, retain_old_id=args.retain_old_id)


