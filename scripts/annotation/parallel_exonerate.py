#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.Annotation import Exonerate

from Routines import AnnotationsRoutines

from CustomCollections.GeneralCollections import IdSet

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with sequences")
parser.add_argument("-a", "--target", action="store", dest="target", required=True,
                    help="File with target sequences")
parser.add_argument("-x", "--annotation", action="store", dest="annotation",
                    help="File with query annotation. Usually cds coordinates in transcript for cdna2genome model")
parser.add_argument("-o", "--output", action="store", dest="output", #required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")
parser.add_argument("-m", "--model", action="store", dest="model", required=True,
                    help="Model to run")
parser.add_argument("-n", "--number_of_results_to_report", action="store",
                    dest="num_of_results_to_report", default=1, type=int,
                    help="Number of results to report for each input sequence")
parser.add_argument("-u", "--num_of_seq_per_file", action="store", dest="num_of_seq_per_file",
                    type=int, default=None,
                    help="Number of sequences per splited input")
parser.add_argument("-d", "--exonerate_dir", action="store", dest="exonerate_dir", default="",
                    help="Directory with exonerate binary")
parser.add_argument("-e", "--num_of_splited_files", action="store", dest="num_of_splited_files",
                    type=int, default=None,
                    help="Number of splited files")
"""
parser.add_argument("-u", "--num_in_seq_per_file", action="store", dest="num_in_seq_per_file",
                    type=int, default=1000,
                    help="Number of sequences per splited input")
"""
"""
parser.add_argument("-r", "--strand", action="store", dest="strand", default="both",
                    help="Strand to consider. Possible variants: both, forward, backward."
                         "Default: both")

parser.add_argument("-e", "--other_options", action="store", dest="other_options",
                    help="Other augustus options")
parser.add_argument("-c", "--augustus_config_dir", action="store", dest="config_dir",
                    help="Augustus config dir")
parser.add_argument("-p", "--pfam_hmm3", action="store", dest="pfam_db",
                    help="Pfam database in hmm3 format")
parser.add_argument("-w", "--swissprot_blast_db", action="store", dest="swissprot_db",
                    help="Blast database of swissprot")
parser.add_argument("-m", "--masking", action="store", dest="masking",
                    help="Gff of bed file with masking of repeats")
"""
args = parser.parse_args()

if args.num_of_seq_per_file and args.num_of_splited_files:
    raise ValueError("Options -u/--num_of_seq_per_file and -e/--num_of_splited_files can't be set simultaneously")

if (not args.num_of_seq_per_file) and (not args.num_of_splited_files):
    args.num_of_splited_files = 10 * args.threads

Exonerate.threads = args.threads
Exonerate.path = args.exonerate_dir
Exonerate.parallel_alignment(args.input, args.target, args.model, num_of_files=args.num_of_splited_files,
                             num_of_recs_per_file=args.num_of_seq_per_file,
                             show_alignment=True, show_sugar=None, show_cigar=None,
                             show_vulgar=None, show_query_gff=True, show_target_gff=True,
                             store_intermediate_files=True,
                             annotation_file=args.annotation,
                             splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                             number_of_results_to_report=args.num_of_results_to_report,
                             converted_output_dir="converted_output")
