#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--protein_alignment", action="store", dest="pep_alignment", required=True,
                    help="File with protein alignment")
parser.add_argument("-c", "--cds_seqs", action="store", dest="cds_seqs", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with cds sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write codon alignment")
parser.add_argument("-a", "--accordance_file", action="store", dest="accordance_file",
                    help="File with CDS to protein id accordance")
parser.add_argument("-f", "--alignment_format", action="store", dest="alignment_format", default="fasta",
                    help="Format of alignments. Default: fasta")
parser.add_argument("-n", "--cds_seqs_format", action="store", dest="cds_format", default="fasta",
                    help="Format of cds sequences. Default: fasta")
parser.add_argument("-i", "--cds_index_file", action="store", dest="cds_index",
                    help="Biopython index of cds files. Default - construct new")
parser.add_argument("-r", "--retain_cds_index", action="store_true", dest="retain_cds_index",
                    help="Retain constructed index after analysis. Default - False")
args = parser.parse_args()

MultipleAlignmentRoutines.get_codon_alignment_from_files(args.pep_alignment, args.cds_seqs, args.output,
                                                         cds2protein_accordance_file=args.accordance_file,
                                                         alignment_format=args.alignment_format,
                                                         nucleotide_sequence_format=args.cds_format,
                                                         cds_index_file=args.cds_index,
                                                         retain_cds_index=args.retain_cds_index)
