#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines
from RouToolPa.Collections.General import IdList


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-g", "--gene_ids", action="store", dest="gene_ids",
                    help="Comma-separated list of IDs for genes to be extracted")
parser.add_argument("-f", "--gene_id_file", action="store", dest="gene_id_file",
                    help="File with ids of genes to be extracted. Ignored - ig ")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output .gff file")


args = parser.parse_args()

value_list = IdList(filename=args.value_file)
AnnotationsRoutines.extract_gff_records_by_description_value(args.input_gff,
                                                             args.output_gff,
                                                             args.field_id_list,
                                                             value_list,
                                                             retain_comments=False)
