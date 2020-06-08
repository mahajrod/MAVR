#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Population import PSMC

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_t", action="store", dest="bam", required=True,
                    help="Bam files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta file with reference")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use. Default: 1 ")

parser.add_argument("-x", "--max_coverage", action="store", dest="max_coverage", type=int,
                    help="Maximum coverage threshold for variant")
parser.add_argument("-i", "--min_coverage", action="store", dest="min_coverage", type=int,
                    help="Minimum coverage threshold for variant")
parser.add_argument("-q", "--min_base_quality", action="store", dest="min_base_quality", default=30, type=int,
                    help="Minimum base quality. Default: 30")
parser.add_argument("-m", "--min_mapping_quality", action="store", dest="min_mapping_quality", default=30, type=int,
                    help="Minimum mapping quality. Default: 30 ")
parser.add_argument("-a", "--mapping_quality_penalty", action="store", dest="mapping_quality_penalty",
                    type=int,
                    help="Penalty for mapping quality of reads having long mismatches. "
                         "Use 50 for BWA. Default: not set")
parser.add_argument("-s", "--min_rms_mapq", action="store", dest="min_rms_mapq", default=30, type=int,
                    help="Minimum RMS mapping quality. Default: 30 ")
args = parser.parse_args()

PSMC.threads = args.threads
PSMC.generate_diploid_fastq_from_bam(args.reference_fasta, args.bam, args.output_prefix,
                                     args.min_coverage, args.max_coverage, split_dir="split/",
                                     min_base_quality=args.min_base_quality,
                                     min_mapping_quality=args.min_mapping_quality,
                                     adjust_mapping_quality=args.mapping_quality_penalty, min_rms_mapq=args.min_rms_mapq)

