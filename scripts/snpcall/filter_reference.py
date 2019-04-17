#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList
from Pipelines import SNPCallPipeline

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference")
parser.add_argument("-l", "--reference_len", action="store", dest="reference_len",
                    help="File with length of scaffolds from reference. Will be calculated if not set")
parser.add_argument("-m", "--masking_gffs", action="store", dest="masking_gffs", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of GFFs with masking")
parser.add_argument("-a", "--annotation_gff", action="store", dest="annotation_gff",
                    help="Gff with annotation")
parser.add_argument("-w", "--white_ids_file", action="store", dest="white_ids_file",
                    help="File with scaffold ids from white list")
parser.add_argument("-b", "--black_ids_file", action="store", dest="black_ids_file",
                    help="File with scaffold ids from black list")
parser.add_argument("-s", "--max_masked_fraction", action="store", dest="max_masked_fraction", default=0.8, type=float,
                    help="Maximum masked fraction of scaffold to retain. Default: 0.8")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-x", "--max_len", action="store", dest="max_len", default=None, type=int,
                    help="Maximum length of scaffold that can be removed. Default: not set")

args = parser.parse_args()

SNPCallPipeline.filter_reference(args.reference, args.masking_gffs, args.output_prefix,
                                 reference_len_file=args.reference_len,
                                 annotation_gff=args.annotation_gff,
                                 max_masked_fraction=args.max_masked_fraction,
                                 white_scaffold_list=IdList(filename=args.white_ids_file) if args.white_ids_file else [],
                                 black_scaffold_list=IdList(filename=args.black_ids_file) if args.black_ids_file else [],
                                 max_length=args.max_len)
