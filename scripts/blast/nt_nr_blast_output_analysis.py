#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Bio import SearchIO
from Bio.SearchIO import QueryResult
from Bio import Entrez
from Routines.File import read_ids
from Routines.Alignment import get_db_ids

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with BLAST results. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with filtered BLAST results. Default: stdout")
parser.add_argument("-f", "--format", action="store", dest="format", default="blast-tab",
                    help="Format of BLAST file. Allowed: blast-tab, blast-xml, blat-psl, hmmer3-tab, hmmer3-domtab. Default: blast-tab")

parser.add_argument("-e", "--max_e_value", action="store", dest="max_e_value", type=float, default=0.001,
                    help="Maximum threshold of e-value. Default: 0.001")
parser.add_argument("-l", "--min_alignment_length", action="store", dest="min_alignment_length", type=int, default=300,
                    help="Minimum threshold  of alignment length. Default: 300")
parser.add_argument("-a", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")
parser.add_argument("-d", "--database", action="store", dest="db", default="nucleotide",
                    help="database to search")
# Not implemented yet
"""
parser.add_argument("-p", "--min_identity_percent", action="store", dest="min_identity_percent", type=float, default=80,
                    help="Minimum threshold of identity percent. Default: 80 %")
parser.add_argument("-b", "--min_bit_score", action="store", dest="min_bit_score", type=float,
                    help="Minimum threshold of bit score. Default: not set")
parser.add_argument("-m", "--max_mismatch_number", action="store", dest="max_mismatch_number", type=int,
                    help="Maximum number of mismatches. Default: not set")
parser.add_argument("-g", "--max_gap_number", action="store", dest="max_gap_number", type=int,
                    help="Maximum number of gaps. Default: not set")
"""
args = parser.parse_args()

Entrez.email = args.email

allowed_formats = ["blast-tab",
                   "blast-xml",
                   "blat-psl",
                   "hmmer3-tab",
                   "hmmer3-domtab"]

if args.format not in allowed_formats:
    raise ValueError("This format of input file is not allowed")

"""
def expression_hsp(hsp):
    # hit_span - length of hit for single-fragment HSP(blast etc). DO NOT work with exonerate
    return (hsp.evalue <= args.max_e_value) and (hsp.hit_span >= args.min_alignment_length)


def iterator(blast_results):
    for query_id in blast_results:
        filtered_query = blast_results[query_id].hsp_filter(func=expression_hsp)
        if filtered_query:
            yield filtered_query
"""

blast_results = SearchIO.index(args.input, args.format)

print(get_db_ids(blast_results))
#SearchIO.write(iterator(blast_results), args.output, args.format)
