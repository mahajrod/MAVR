#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import SNPeff

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input table")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output table with common name")
parser.add_argument("-s", "--synonym_file", action="store", dest="synonym_file", required=True,
                    help="File with synonyms")
parser.add_argument("-e", "--header_name_for_synonym", action="store", dest="header_name_for_synonym",
                    default="Common_name",
                    help="Header name for synonym column")
parser.add_argument("-k", "--key_index", action="store", dest="key_index", type=int, default=0,
                    help="Key column in file with synonyms(0-based). Default: 0")
parser.add_argument("-v", "--value_index", action="store", dest="value_index", type=int, default=1,
                    help="Value column in file with synonyms(0-based). Default: 1")
parser.add_argument("-a", "--snpeff_tab_column_id_column", action="store", dest="snpeff_tab_column_id_column", type=int, default=8,
                    help="SNPeff tab file ID column. Default: 0")

args = parser.parse_args()

SNPeff.add_gene_synonyms(args.input, args.output, args.synonym_file,
                         key_column=args.key_index, value_column=args.value_index,
                         header_name_for_synonym=args.header_name_for_synonym,
                         snpeff_tab_column_id_column=args.snpeff_tab_column_id_column)
