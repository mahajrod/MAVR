#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.BLAST import DustMasker, MakeBLASTDb, BLASTn, BLASTp

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-m", "--mask", action="store_true", dest="mask",
                    help="Mask sequences with dustmasker")
parser.add_argument("-n", "--name", action="store", dest="name", default="tmp_database",
                    help="Database name")
parser.add_argument("-y", "--database_type", action="store", dest="database_type", default="nucleotide",
                    help="Database type. Allowed types: protein and nucleotide. Default: nucleotide")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use")
parser.add_argument("-b", "--other_blast_options", action="store", dest="other_options",
                    help="Other blast options")
parser.add_argument("-e", "--evalue", action="store", dest="evalue",
                    help="E-value cutoff")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="E-value cutoff")
parser.add_argument("-f", "--output_format", action="store", dest="output_format", default=6, type=int,
                    help="Output format:"
                         "0 = pairwise,"
                         "1 = query-anchored showing identities,"
                         "2 = query-anchored no identities,"
                         "3 = flat query-anchored, show identities,"
                         "4 = flat query-anchored, no identities,"
                         "5 = XML Blast output,"
                         "6 = tabular,"
                         "7 = tabular with comment lines,"
                         "8 = Text ASN.1,"
                         "9 = Binary ASN.1,"
                         "10 = Comma-separated values,"
                         "11 = BLAST archive format (ASN.1)")
args = parser.parse_args()

if args.mask:
    mask_file = args.input + ".asnb"
    DustMasker.mask(args.input_file, mask_file)

if args.database_type == "nucleotide":
    MakeBLASTDb.make_nucleotide_db(args.input, args.name, mask_file if args.mask else None,
                                   output_file=args.name)
    BLASTn.parallel_blastn(args.input, args.name, outfile=args.output,
                           blast_options=args.other_options, split_dir="splited_fasta",
                           splited_output_dir="splited_output_dir",
                           evalue=args.evalue, output_format=args.output_format,
                           threads=args.threads,
                           combine_output_to_single_file=True)
elif args.database_type == "protein":
    MakeBLASTDb.make_protein_db(args.input, args.name, mask_file if args.mask else None,
                                output_file=args.name)
    BLASTp.parallel_blastp(args.input, args.name, outfile=args.output,
                           blast_options=args.other_options, split_dir="splited_fasta",
                           splited_output_dir="splited_output_dir",
                           evalue=args.evalue, output_format=args.output_format,
                           threads=args.threads,
                           combine_output_to_single_file=True)

