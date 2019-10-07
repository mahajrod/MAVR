#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Pipelines import AlignmentPipeline

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_dir", action="store", dest="sample_dir", required=True,
                    help="Directory with samples data")
parser.add_argument("-o", "--outdir", action="store", dest="outdir", default="./",
                    help="Output directory. Default: current directory")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",  type=lambda s: s.split(","),
                    help="List of samples to call variants for. Default: all samplews in sample directory")
parser.add_argument("-i", "--aligner_index", action="store", dest="index",
                    help="Aligner-specific index")
parser.add_argument("-a", "--aligner", action="store", dest="aligner", default="bwa",
                    help="Aligner to use. Possible aligners: bwa(default), bowtie2")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-q", "--quality", action="store", dest="quality", default="phred33",
                    help="Quality type. Possible variants phred33, phred64. Default: phred33")
parser.add_argument("-j", "--alignment_format", action="store", dest="alignment_format", default="bam",
                    help="Format of output alignments. Allowed: bam(default), sam, cram")
parser.add_argument("-g", "--add_read_groups_by_picard", action="store_true",
                    dest="add_read_groups_by_picard", default=False,
                    help="Add read groups to final bam using PICARD. Use this option "
                         "if aligner don't support adding readgroups itself")
parser.add_argument("-c", "--picard_dir", action="store", dest="picard_dir", default="",
                    help="Path to Picard directory. Required to add read groups and mark duplicates")
parser.add_argument("-l", "--aligner_dir", action="store", dest="aligner_dir", default="",
                    help="Path to aligner directory")
parser.add_argument("-u", "--read_file_suffix", action="store", dest="read_file_suffix", default="",
                    help="Suffix of read files, i.e forward read file have folling name"
                         " <sample><read_suffix>_1.<extension> Default: ''")
parser.add_argument("-x", "--read_file_extension", action="store", dest="read_file_extension", default="fastq",
                    help="Extension of read files, i.e forward read file have folling name"
                         " <sample><read_suffix>_1.<extension> Default: 'fastq'")
parser.add_argument("-z", "--gzipped_reads", action="store_true", dest="gzipped_reads", default=False,
                    help="Reads are gzipped")
parser.add_argument("-k", "--skip_duplicates", action="store_true", dest="skip_duplicates", default=False,
                    help="Skip duplicate marking")
"""
parser.add_argument("-z", "--calculate_median_coverage", action="store_true", dest="calculate_median_coverage",
                    default=False,
                    help="Calculate median coverage")
parser.add_argument("-x", "--calculate_mean_coverage", action="store_true", dest="calculate_mean_coverage",
                    default=False,
                    help="Calculate mean coverage")
parser.add_argument("-f", "--flanks_size", action="store", dest="flanks_size",
                    default=0, type=int,
                    help="Size of flanks to remove when calculating mean/median coverage")
"""


parser.add_argument("-m", "--max_memory", action="store", dest="max_memory", default=30, type=int,
                    help="Maximum memory to use in gigabytes. Default: 30")

args = parser.parse_args()

AlignmentPipeline.max_memory = args.max_memory
AlignmentPipeline.threads = args.threads
AlignmentPipeline.BWA_dir = args.aligner_dir
AlignmentPipeline.bowtie2_dir = args.aligner_dir
AlignmentPipeline.Picard_dir = args.picard_dir

AlignmentPipeline.align(args.sample_dir, args.index, aligner=args.aligner, sample_list=args.sample_list,
                        outdir=args.outdir, quality_score_type=args.quality, read_suffix=args.read_file_suffix,
                        read_extension=args.read_file_extension, alignment_format=args.alignment_format,
                        threads=None, mark_duplicates=not args.skip_duplicates, platform="Illumina",
                        add_read_groups_by_picard=args.add_read_groups_by_picard, gzipped_reads=args.gzipped_reads)
