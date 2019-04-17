#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cds_file", action="store", dest="cds_file", required=True,
                    help="Input file CDS sequences")
parser.add_argument("-p", "--pep_file", action="store", dest="pep_file", required=True,
                    help="Input file protein sequences")
parser.add_argument("-o", "--output_file", action="store", dest="out", required=True,
                    help="Output file")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed: fasta, genbank. Default: fasta")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print warning if no protein was found for CDS")
parser.add_argument("-m", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode of sequence files. Allowed: parse, index, index_db."
                         "Default: parse")
parser.add_argument("-t", "--genetic_code_table", action="store", dest="genetic_code_table", default=1, type=int,
                    help="Genetic code to use for translation of CDS. "
                         "Allowed: table number from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
                         "Default: 1(The standard code)")
parser.add_argument("-d", "--id_check", action="store_true", dest="id_check",
                    help="Also use id check - if there is id present in both files consider them as accordance")
parser.add_argument("-w", "-cds_with_no_pep_idfile", action="store", dest="cds_with_no_pep_idfile",
                    help="File to write ids of CDS with no protein hit. Default: not set")

args = parser.parse_args()

SequenceRoutines.get_cds_to_pep_accordance_from_files(args.cds_file, args.pep_file, args.out, verbose=args.verbose,
                                                      parsing_mode=args.parsing_mode,
                                                      genetic_code_table=args.genetic_code_table,
                                                      include_id_check=args.id_check,
                                                      cds_with_no_pep_idfile=args.cds_with_no_pep_idfile)
if args.parsing_mode == "index_db":
    os.remove("cds_tmp.idx")
    os.remove("pep_tmp.idx")


