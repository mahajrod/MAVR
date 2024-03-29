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
                    help="List of samples to call variants for. Default: all samples in sample directory")
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
parser.add_argument("-e", "--mkdup_tool", action="store", dest="mkdup_tool", default="samtools",
                    help="Tool to use for duplicate marking. Allowed: samtools(default), sambamba, picard")
parser.add_argument("-c", "--picard_dir", action="store", dest="picard_dir", default="",
                    help="Path to Picard directory. Required to add read groups and mark duplicates")
parser.add_argument("-b", "--sambamba_dir", action="store", dest="sambamba_dir", default="",
                    help="Path to sambamba directory. Required to add mark duplicates")
parser.add_argument("-l", "--aligner_dir", action="store", dest="aligner_dir", default="",
                    help="Path to aligner directory")
parser.add_argument("--local_aln", action="store_true", dest="local_aln", default=False,
                    help="Perform local alignment. Affects only Bowtie2 as default mode for Bowtie2 is end-to-end."
                         "Default: False")
parser.add_argument("-u", "--read_file_suffix", action="store", dest="read_file_suffix", default="",
                    help="Suffix of read files, i.e forward read file have folling name"
                         " <sample><read_suffix>_1.<extension> Default: ''")
parser.add_argument("--unpaired_read_file_suffix", action="store", dest="unpaired_read_file_suffix", default=None,
                    help="Suffix of unpaired read files, i.e forward unpaired file have folling name"
                         " <sample><read_suffix> Default: None")
parser.add_argument("-x", "--read_file_extension", action="store", dest="read_file_extension", default="fastq",
                    help="Extension of read files, i.e forward read file have folling name"
                         " <sample><read_suffix>_1.<extension> Default: 'fastq'")
parser.add_argument("-z", "--gzipped_reads", action="store_true", dest="gzipped_reads", default=False,
                    help="Reads are gzipped")
parser.add_argument("-k", "--skip_duplicates", action="store_true", dest="skip_duplicates", default=False,
                    help="Skip duplicate marking")
parser.add_argument("-r", "--retain_intermediate_files", action="store_true", dest="retain_intermediate_files",
                    default=False,
                    help="Retain intermediate files. Default: False")
parser.add_argument("--tmp_dir", action="store", dest="tmp_dir",
                    help="Directory for temporary files. Default: system default")
parser.add_argument("-w", "--calculate_coverage", action="store_true", dest="calculate_coverage", default=False,
                    help="Calculate coverage. Default: False")
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
parser.add_argument("--softclipping_penalty", action="store", dest="softclipping_penalty",
                    type=lambda s: map(int, s.split(",")),
                    help="Softclipping penalty for 5' and 3' end, separated by comma, i.e 5,5. Default: bwa default")
args = parser.parse_args()

AlignmentPipeline.max_memory = args.max_memory
AlignmentPipeline.threads = args.threads
AlignmentPipeline.BWA_dir = args.aligner_dir
AlignmentPipeline.bowtie2_dir = args.aligner_dir
AlignmentPipeline.Picard_dir = args.picard_dir
AlignmentPipeline.sambamba_dir = args.sambamba_dir
AlignmentPipeline.tmp_dir = args.tmp_dir


AlignmentPipeline.align(args.sample_dir, args.index, aligner=args.aligner, sample_list=args.sample_list,
                        outdir=args.outdir, quality_score_type=args.quality, read_suffix=args.read_file_suffix,
                        read_extension=args.read_file_extension, alignment_format=args.alignment_format,
                        unpaired_read_suffix=args.unpaired_read_file_suffix,
                        threads=None, mark_duplicates=not args.skip_duplicates, platform="Illumina",
                        add_read_groups_by_picard=args.add_read_groups_by_picard, gzipped_reads=args.gzipped_reads,
                        keep_inremediate_files=args.retain_intermediate_files,
                        mark_duplicates_tool=args.mkdup_tool,
                        calculate_coverage=args.calculate_coverage, softclipping_penalty=args.softclipping_penalty,
                        local_alignment=args.local_aln)
