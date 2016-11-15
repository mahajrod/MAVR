#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Pipelines import DiffExpressionPipeline

from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-g", "--genome_dir", action="store", dest="genome_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with star index for genome")
parser.add_argument("-f", "--genome_fasta", action="store", dest="genome_fasta",
                    type=os.path.abspath,
                    help="Path to genome fasta file. If set Star will construct genome index first"
                         "in directory set by -g/--genome_dir")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use in Trimmomatic. Default - 1.")
parser.add_argument("-a", "--annotation_gtf", action="store", dest="annotation_gtf", type=os.path.abspath,
                    help="Gtf file with annotations")
parser.add_argument("-i", "--genome_size", action="store", dest="genome_size", type=int,
                    help="Genome size. Required for constructing genome index")
parser.add_argument("-j", "--junction_tab_file", action="store", dest="junction_tab_file",
                    help="Junction tab file")
parser.add_argument("-r", "--star_dir", action="store", dest="star_dir", default="",
                    help="Directory with STAR binary")
parser.add_argument("-u", "--include_unmapped_reads", action="store_true",
                    dest="include_unmapped_reads", default=False,
                    help="Include unmapped reads in Bam file")


""""
parser.add_argument("-x", "--general_stat_file", action="store", dest="general_stat_file", required=True,
                    help="File to write general statistics about filtration")
parser.add_argument("-k", "--adapter_kmers", action="store", dest="adapter_kmers", type=os.path.abspath,
                    required=True,
                    help="File with adapter k-mers for Coockiecutter")

parser.add_argument("-m", "--mismatch_number", action="store", dest="mismatch_number", type=int, default=2,
                    help="Number of mismatches in adapter seed. Works only if -a/--adapters option is set. Default - 2.")

parser.add_argument("-p", "--pe_score", action="store", dest="pe_score", type=int, default=30,
                    help="PE reads adapter score. Works only if -a/--adapters option is set. Default - 30.")
parser.add_argument("-e", "--se_score", action="store", dest="se_score", type=int, default=10,
                    help="SE reads adapter score. Works only if -a/--adapters option is set. Default - 10.")
parser.add_argument("-n", "--min_adapter_len", action="store", dest="min_adapter_len", type=int, default=1,
                    help="Minimum length of adapter fragment. Works only if -a/--adapters option is set. Default - 1.")

parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=20,
                    type=int,
                    help="Quality threshold for sliding window or whole read."
                         "Depends on -q/--average_quality_threshold option.Default - 15.")
parser.add_argument("-b", "--base_quality", action="store", dest="base_quality", default="phred33",
                    help="Type of base quality. Possible variants: phred33, phred64. Default - phred33 ")

parser.add_argument("-l", "--min_length", action="store", dest="min_len", type=int, default=50,
                    help="Minimum length of read to retain. Default - 50")


parser.add_argument("-c", "--coockiecutter_dir", action="store", dest="coockiecutter_dir", default="",
                    help="Path to Coockiecutter directory")

"""

args = parser.parse_args()

"""
EXAMPLE
skliver@supermicro:
cd ~/workdir/yeast/nizhnikov/good_run/

"""

DiffExpressionPipeline.star_and_htseq(args.genome_dir, args.samples_dir, args.output_dir, genome_fasta=args.genome_fasta,
                                      genome_size=args.genome_size, samples_to_handle=args.samples,
                                      annotation_gtf=args.annotation_gtf,
                                      feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                                      exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                                      junction_tab_file_list=args.junction_tab_file,
                                      three_prime_trim=None, five_prime_trim=None,
                                      adapter_seq_for_three_prime_clip=None,
                                      max_mismatch_percent_for_adapter_trimming=None,
                                      three_prime_trim_after_adapter_clip=None,
                                      output_type="BAM", sort_bam=True, max_memory_for_bam_sorting=None,
                                      include_unmapped_reads_in_bam=args.include_unmapped_reads,
                                      output_unmapped_reads=True,  two_pass_mode=False, star_dir=args.star_dir,
                                      threads=args.threads)
