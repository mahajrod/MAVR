#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input gff from AUGUSTUS")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff with replaced ids")
parser.add_argument("-g", "--gene_syn_file", action="store", dest="gene_syn_file", required=True,
                    help="File with gene synonyms")
parser.add_argument("-t", "--transcript_syn_file", action="store", dest="transcript_syn_file", required=True,
                    help="File with transcript synonyms")
parser.add_argument("-c", "--cds_syn_file", action="store", dest="cds_syn_file",
                    help="File with CDS synonyms")

args = parser.parse_args()

AUGUSTUS.replace_augustus_ids_by_syn(args.input_gff, args.output_gff, args.gene_syn_file,
                                     args.transcript_syn_file, cds_syn_file=args.cds_syn_file)
