#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Pipelines import ITSPipeline
from RouToolPa.Routines.File import check_path


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-s", "--samples", action="store", dest="samples", type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-f", "--filtered", action="store_true", dest="filtered", default=False,
                    help="Input reads were already filtered and filenames follow convention."
                         "Use this option if you have previously ran pipeline and wish to try new reference"
                         "Default: False")
parser.add_argument("--reference", action="store", dest="reference", required=True,
                    help="Reference fasta file")
parser.add_argument("--index", action="store", dest="index", required=True,  help="BWA index")

parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("--output_prefix", action="store", dest="output_prefix", required=True,
                    default="./", help="Prefix for output vcf file")

parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use in Trimmomatic. Default - 1.")

parser.add_argument("-a", "--adapters", action="store", dest="adapters", type=os.path.abspath,
                    help="File with adapters to trim by Trimmomatic")
parser.add_argument("-k", "--adapter_kmers", action="store", dest="adapter_kmers", type=os.path.abspath,
                    help="File with adapter k-mers ")

parser.add_argument("-m", "--mismatch_number", action="store", dest="mismatch_number", type=int, default=2,
                    help="Number of mismatches in adapter seed. Works only if -a/--adapters option is set. Default - 2.")
parser.add_argument("-p", "--pe_score", action="store", dest="pe_score", type=int, default=30,
                    help="PE reads adapter score. Works only if -a/--adapters option is set. Default - 30.")
parser.add_argument("-e", "--se_score", action="store", dest="se_score", type=int, default=10,
                    help="SE reads adapter score. Works only if -a/--adapters option is set. Default - 10.")
parser.add_argument("-n", "--min_adapter_len", action="store", dest="min_adapter_len", type=int, default=1,
                    help="Minimum length of adapter fragment. Works only if -a/--adapters option is set. Default - 1.")
parser.add_argument("-g", "--sliding_window_size", action="store", dest="sliding_window_size", type=int,
                    help="Size of sliding window when checking quality. "
                         "If not set - reads will be filtered by mean quality")
parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=20,
                    type=int,
                    help="Quality threshold for sliding window or whole read."
                         "Depends on -q/--average_quality_threshold option.Default - 20.")
parser.add_argument("-b", "--base_quality", action="store", dest="base_quality", default="phred33",
                    help="Type of base quality. Possible variants: phred33, phred64. Default - phred33 ")
parser.add_argument("-l", "--min_length", action="store", dest="min_len", type=int, default=50,
                    help="Minimum length of read to retain. Default - 50")
parser.add_argument("-j", "--trimmomatic_dir", action="store", dest="trimmomatic_dir", default="",
                    help="Path to Trimmomatic directory")
parser.add_argument("-c", "--trimmer_dir", action="store", dest="trimmer_dir", default="",
                    help="Path to Trimmer directory")
parser.add_argument("-w", "--bamutil_dir", action="store", dest="bamutil_dir", default="",
                    help="Path to BamUtil directory")
parser.add_argument("-r", "--keep_intermediate_files", action="store_true",
                    dest="keep_intermediate_files", default=False,
                    help="Keep intermediate files.Default: False")

args = parser.parse_args()

"""
EXAMPLE
~/Soft/MAVR/scripts/ITS/ITS_pipeline.py -d ../reads/raw/ \
                                        -s A01,A02 \
                                        --reference ../reference/acipenser_ruthenus.ITS.fasta \
                                        --index ../reference/acipenser_ruthenus.ITS.fasta \
                                        -o ./ -t 20 \
                                        -a ~/Soft/Trimmomatic-0.36/adapters/TruSeq2-PE.fa \
                                        -k ~/data/service_seq/illumina_adapters_with_rev_com_23_mer.kmer \
                                        --output_prefix acipenser_ruthenus.ITS \
                                        -j ~/Soft/Trimmomatic-0.36/
                                        
~/Soft/MAVR/scripts/ITS/ITS_pipeline.py -f -d ../acipenser_ruthenus/reads/filtered/ \
                                        --reference ../reference/acipenser_baerii.ITS.fasta \
                                        --index ../reference/acipenser_baerii.ITS.fasta \
                                        -o ./ -t 20 \
                                        -a ~/Soft/Trimmomatic-0.36/adapters/TruSeq2-PE.fa \
                                        -k ~/data/service_seq/illumina_adapters_with_rev_com_23_mer.kmer \
                                        --output_prefix acipenser_ruthenus.ITS \
                                        -j ~/Soft/Trimmomatic-0.36/
                                        
~/Soft/MAVR/scripts/ITS/ITS_pipeline.py -d ../reads/raw/ \
                                        --reference ../reference/acipenser_ruthenus.ITS.fasta \
                                        --index ../reference/acipenser_ruthenus.ITS.fasta \
                                        -o ./ -t 20 \
                                        -a ~/Soft/Trimmomatic-0.36/adapters/TruSeq2-PE.fa \
                                        -k ~/data/service_seq/illumina_adapters_with_rev_com_23_mer.kmer \
                                        --output_prefix acipenser_ruthenus.ITS \
                                        -j ~/Soft/Trimmomatic-0.36/                                       
"""

ITSPipeline.pipeline(args.samples_dir, args.output_dir, args.adapter_kmers, args.adapters,
                     args.reference, args.index, args.output_prefix,
                     samples_to_handle=args.samples, threads=args.threads, trimmomatic_dir=args.trimmomatic_dir,
                     trimmer_dir=args.trimmer_dir, bam_util_dir=args.bamutil_dir,
                     mismatch_number=args.mismatch_number, pe_reads_score=args.pe_score,
                     se_read_score=args.se_score, min_adapter_len=args.min_adapter_len,
                     sliding_window_size=args.sliding_window_size,
                     average_quality_threshold=args.average_quality_threshold,
                     leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
                     crop_length=None, head_crop_length=None, min_len=args.min_len,
                     base_quality=args.base_quality,
                     remove_intermediate_files=not args.keep_intermediate_files, filtered_reads=args.filtered)


