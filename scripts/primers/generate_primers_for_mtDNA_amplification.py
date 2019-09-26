#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Pipelines import MitochondrialAmplificationPrimerPipeline

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mitochondrial_fasta", action="store", dest="mitochondrial_fasta", required=True,
                    help="Fasta file with mtDNA sequences")
parser.add_argument("-g", "--genome_fasta", action="store", dest="genome_fasta", required=True,
                    help="Fasta file with whole genome sequence")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-c", "--coordinates", action="store", dest="coordinates",
                    help="File with coordinates of regions to amplify")
parser.add_argument("--max_pcr_product_len", action="store", dest="max_pcr_product_len", type=int, default=4700,
                    help="Maximum length of PCR product")
parser.add_argument("--min_pcr_product_len", action="store", dest="min_pcr_product_len", type=int, default=4200,
                    help="Minimumlength of PCR product")

parser.add_argument("-k", "--directory_with_kmer_counts", action="store", dest="directory_with_kmer_counts", required=True,
                    help="Directory with files containing kmer counts")
parser.add_argument("-x", "--kmer_file_prefix", action="store", dest="kmer_file_prefix", required=True,
                    help="Prefix of files with kmer counts")
parser.add_argument("-z", "--count_kmers", action="store_true", dest="count_kmers",
                    help="Count kmers by Glistmaker for primer 3 using path specified "
                         "by -k/--directory_with_kmer_counts and -x/--kmer_file_prefix to store the results")


parser.add_argument("--primer3_dir", action="store", dest="primer3_dir", default="",
                    help="Directory with primer3_core binary")
parser.add_argument("--glistmaker_dir", action="store", dest="glistmaker_dir", default="",
                    help="Directory with Glistmaker binary")
parser.add_argument("--trf_dir", action="store", dest="trf_dir", default="",
                    help="Directory with TRF binary")

parser.add_argument("--primer3_thermo_config_dir", action="store", dest="primer3_thermo_config_dir",
                    help="Directory with primer3 config for thermodynamic approach")

parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default: 1")
"""
parser.add_argument("-s", "--split_output_by_monomer_len", action="store_true",
                    dest="split_output_by_monomer_len", default=False,
                    help="Split output by STR monomer length")
"""
args = parser.parse_args()

MitochondrialAmplificationPrimerPipeline.primer_prediction_pipeline(args.coordinates,
                                                                    args.mitochondrial_fasta,
                                                                    args.genome_fasta,
                                                                    args.output_prefix,
                                                                    kmer_dir=args.directory_with_kmer_counts,
                                                                    kmer_file_prefix=args.kmer_file_prefix,
                                                                    count_kmers=args.count_kmers,
                                                                    pcr_product_size_range=(args.min_pcr_product_len,
                                                                                            args.max_pcr_product_len),
                                                                    optimal_primer_len=None,
                                                                    min_primer_len=20,
                                                                    max_primer_len=28,
                                                                    max_ns_accepted=None,
                                                                    softmasked_input=False,
                                                                    optimal_GC=None, min_GC=None, max_GC=None,
                                                                    optimal_melting_temperature=69,
                                                                    min_melting_temperature=68,
                                                                    max_melting_temperature=70,
                                                                    black_list_of_seqs_fasta=None,
                                                                    threads=None,)

