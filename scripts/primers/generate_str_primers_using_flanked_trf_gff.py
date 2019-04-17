#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from RouToolPa.Routines import FileRoutines
from Pipelines import STRPrimerPipeline

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--trf_flank_gff", action="store", dest="trf_flank_gff", required=True,
                    help="TRF Gff file with flanks added")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-s", "--fasta_with_flanks", action="store", dest="fasta_with_flanks", required=True,
                    help="Fasta file with flanked repeat sequences")
parser.add_argument("-k", "--directory_with_kmer_counts", action="store", dest="directory_with_kmer_counts", required=True,
                    help="Directory with files containing kmer counts")
parser.add_argument("-r", "--kmer_file_prefix", action="store", dest="kmer_file_prefix", required=True,
                    help="Prefix of files with kmer counts")

parser.add_argument("-p", "--primer3_dir", action="store", dest="primer3_dir", default="",
                    help="Directory with primer3_core binary")
parser.add_argument("-y", "--primer3_thermo_config_dir", action="store", dest="primer3_thermo_config_dir",
                    help="Directory with primer3 config for thermodynamic approach")
parser.add_argument("-m", "--format_output", action="store_true", dest="format_output", default=False,
                    help="Convert output to human readable form")

"""
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default: 1")

parser.add_argument("-p", "--glistmaker_path", action="store", dest="glistmaker_path",
                    type=FileRoutines.split_filename, default=["", "glistmaker", ""],
                    help="Path to Glistmaker binary")
"""

args = parser.parse_args()

STRPrimerPipeline.primer3_dir = args.primer3_dir
STRPrimerPipeline.predict_primers(args.trf_flank_gff, args.fasta_with_flanks, args.output_prefix,
                                  args.directory_with_kmer_counts, args.kmer_file_prefix, pcr_product_size_range=None,
                                  optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                                  softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                                  optimal_melting_temperature=None, min_melting_temperature=None,
                                  max_melting_temperature=None, black_list_of_seqs_fasta=None,
                                  thermodynamic_parameters_dir=args.primer3_thermo_config_dir,
                                  format_output=args.format_output)

