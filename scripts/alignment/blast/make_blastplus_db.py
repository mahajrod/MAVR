#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.BLAST import *

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file")
parser.add_argument("-m", "--mask", action="store_true", dest="mask",
                    help="Mask sequences with dustmasker")
parser.add_argument("-t", "--database_title", action="store", dest="db_title",
                    help="Title of database")
parser.add_argument("-n", "--name", action="store", dest="name",
                    help="Database name")
parser.add_argument("-y", "--database_type", action="store", dest="database_type", default="nucleotide",
                    help="Database type. Allowed types: protein and nucleotide. Default: nucleotide")
args = parser.parse_args()

if args.mask:
    mask_file = args.input_file + ".asnb"
    DustMasker.mask(args.input_file, mask_file)

if args.database_type == "nucleotide":
    MakeBLASTDb.make_nucleotide_db(args.input_file, args.db_title, mask_file if args.mask else None,
                                   output_file=args.name)

elif args.database_type == "protein":
    MakeBLASTDb.make_protein_db(args.input_file, args.db_title, mask_file if args.mask else None,
                                output_file=args.name)

BLASTDbCmd.db_info(args.name)



