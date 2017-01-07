#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Pipelines.GenomeAssembly import ScaffoldingPipeline

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward", action="store", dest="forward", required=True,
                    help="Comma separated list of files with forward sequences")
parser.add_argument("-r", "--reverse", action="store", dest="reverse", required=True,
                    help="Comma separated list of files with reverse sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-i", "--estimated_insert_size", action="store", dest="estimated_insert_size", required=True,
                    type=int, help="Estimated insert size")
parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                    help="Path to file with genome")
parser.add_argument("-b", "--bowtie2_index", action="store", dest="bowtie2_index", required=True,
                    help="Path to bowtie2 index of genome")
parser.add_argument("-e", "--read_orientation", action="store", dest="read_orientation", default="fr",
                    help="Read orientation in pair. Allowed: fr(illumina paired end, default), "
                         "rf(illumina mate-pairs), ff")
parser.add_argument("-m", "--format", action="store", dest="format", default="fasta",
                    help="Format of genome file. Allowed formats genbank, fasta(default)")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")

args = parser.parse_args()


"""
example of usage

~/Soft/MAVR/scripts/sequence/check_pairing.py -p parse \
                                              -a ".F" \
                                              -b ".R" \
                                              -o GSS_BOH_BAC_end \
                                              -f GSS_BOH_BAC_end.forward.fa \
                                              -r GSS_BOH_BAC_end.reverse.fa

"""
ScaffoldingPipeline.get_insert_size_distribution(os.getcwd(), args.forward, args.reverse_files,
                                                 args.estimated_insert_size, args.output_prefix,
                                                 args.genome, args.bowtie2_index, read_orientation="fr",
                                                 parsing_mode=args.parsing_mode, number_of_bins=100,
                                                 genome_format=args.format)
