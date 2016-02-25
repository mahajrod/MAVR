#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import AlignmentRoutines, FileRoutines


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

args = parser.parse_args()

AlignmentRoutines.get_codon_alignment_from_files(args.pep_alignment, args.cds_seqs, args.output,
                                                 cds2protein_accordance_file=args.accordance_file,
                                                 alignment_format=args.alignment_format,
                                                 nucleotide_sequence_format=args.cds_format)
