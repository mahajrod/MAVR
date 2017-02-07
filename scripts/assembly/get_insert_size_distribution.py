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
parser.add_argument("-b", "--genome_index", action="store", dest="genome_index", required=True,
                    help="Path to genome index of genome. Must be compatible with choosed aligner")
parser.add_argument("-e", "--read_orientation", action="store", dest="read_orientation", default="fr",
                    help="Read orientation in pair. Allowed: fr(illumina paired end, default), "
                         "rf(illumina mate-pairs), ff. WORKS WITH BOWTIE2 ONLY!!!")
parser.add_argument("-m", "--format", action="store", dest="format", default="fasta",
                    help="Format of genome file. Allowed formats genbank, fasta(default)")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")
parser.add_argument("-a", "--input_fasta", action="store_true", dest="input_fasta",
                    help="Reads are in fasta format. Default: False. NEEDS FOR BOWTIE2 ONLY!!!")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")
parser.add_argument("-d", "--aligner_binary_dir", action="store", dest="aligner_binary_dir", default="",
                    help="Directory with aligner binary")
parser.add_argument("-s", "--store_sam", action="store_true", dest="store_sam",
                    help="Store .sam file with read alignments")
parser.add_argument("-l", "--aligner", action="store", dest="aligner", default="bwa",
                    help="Aligner. Allowed: bowtie2, bwa(default)")
parser.add_argument("-x", "--max_xlimit_for_histo", action="store",
                    dest="xlimit_for_histo",
                    help="Xlimit for histogram. Default: 3 * estimated_insert_size "
                         "(set by -i/--estimated_insert_size option) ")

args = parser.parse_args()


"""
example of usage

cd /mnt/guatemala/skliver/Boechera/Boechera_holboellii/genome/insert_size_estimation/BES$

~/Soft/MAVR/scripts/assembly/get_insert_size_distribution.py -f ../../gz/GSS_BOH_BAC_end.pe.forward.fasta \
                                                             -r ../../gz/GSS_BOH_BAC_end.pe.reverse.fasta \
                                                             -o GSS_BOH_BAC.rf \
                                                             -i 100000 \
                                                             -e rf \
                                                             -p parse \
                                                             -b /mnt/peru/skliver/Boechera/Boechera_holboellii/genome/assemblies/discovar_adapter_filtered_cookie_trimmomatic/a.final/bowtie2_index/a.lines -g /mnt/peru/skliver/Boechera/Boechera_holboellii/genome/assemblies/discovar_adapter_filtered_cookie_trimmomatic/a.final/a.lines.fasta \
                                                             -a \
                                                             -t 30

"""
ScaffoldingPipeline.threads = args.threads
ScaffoldingPipeline.get_insert_size_distribution(os.getcwd(), args.forward, args.reverse,
                                                 args.estimated_insert_size, args.output_prefix,
                                                 args.genome, args.genome_index, read_orientation="fr",
                                                 parsing_mode=args.parsing_mode, number_of_bins=100,
                                                 genome_format=args.format, input_files_are_fasta=args.input_fasta,
                                                 store_sam=args.store_sam, aligner=args.aligner,
                                                 aligner_binary_dir=args.aligner_binary_dir,
                                                 xlimit_for_histo=args.xlimit_for_histo)
