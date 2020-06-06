#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Samtools import VariantCall

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_list", action="store", dest="bam_list", type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of bam files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta file with reference")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use. Default: 1 ")
parser.add_argument("-k", "--chunk_length", action="store", dest="chunk_length", default=1000000, type=int,
                    help="Chunk length. Default: 1 000 000 bp")
parser.add_argument("-d", "--max_per_sample_coverage", action="store", dest="max_per_sample_coverage", type=int,
                    help="Maximum per sample coverage to use")
parser.add_argument("-q", "--min_base_quality", action="store", dest="min_base_quality", default=30, type=int,
                    help="Minimum base quality. Default: 30")
parser.add_argument("-m", "--min_mapping_quality", action="store", dest="min_mapping_quality", default=30, type=int,
                    help="Minimum mapping quality. Default: 30 ")
parser.add_argument("-a", "--mapping_quality_penalty", action="store", dest="mapping_quality_penalty",
                    type=int,
                    help="Penalty for mapping quality of reads having long mismatches. "
                         "Use 50 for BWA. Default: not set")
parser.add_argument("-c", "--consensus_caller_model", action="store_true", dest="consensus_caller_model",
                    help="Use consensus caller model(old one implemented in samtools). Default: false")

args = parser.parse_args()


VariantCall.threads = args.threads
VariantCall.call_variants(args.reference, args.output_prefix, args.bam_list, chunk_length=args.chunk_length,
                          split_dir="split/", max_coverage=args.max_per_sample_coverage,
                          min_base_quality=args.min_base_quality, min_mapping_quality=args.min_mapping_quality,
                          adjust_mapping_quality=args.mapping_quality_penalty,
                          consensus_caller_model=args.consensus_caller_model)
