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
parser.add_argument("--aligned", action="store_true", dest="aligned", default=False,
                    help="Input reads were already aligned and filenames follow convention."
                         "Use this option if you have previously ran pipeline "
                         "and wish to try new variant calling options"
                         "Default: False")
parser.add_argument("--aligned_and_clipped", action="store_true", dest="aligned_and_clipped", default=False,
                    help="Input reads were already aligned, clipped and filenames follow convention."
                         "Use this option if you have previously ran pipeline "
                         "and wish to try new variant calling options"
                         "Default: False")
parser.add_argument("--reference", action="store", dest="reference", required=True,
                    help="Reference fasta file")
parser.add_argument("--max_coverage_for_variant_call", action="store", dest="max_coverage_for_variant_call", default=0,
                    type=int,
                    help="Maximum reads to use for each position. Default: 0, i.e not limited")
parser.add_argument("--min_coverage_for_filtering", action="store", dest="min_coverage_for_filtering", default=100,
                    type=int,
                    help="Minimum coverage to retain variant during filtering. Default: 100")
parser.add_argument("--adjust_mapping_quality", action="store", dest="adjust_mapping_quality", default=None,
                    type=int,
                    help="Adjust mapping quality of reads. Coefficient for downgrading mapping quality for "
                         "reads containing excessive mismatches. Given a read with a phred-scaled probability q of "
                         "being generated from the mapped position, "
                         "the new mapping quality is about sqrt((INT-q)/INT)*INT. "
                         "A zero value disables this functionality; if enabled, the recommended value for BWA is 50. "
                         "Default: not set")
parser.add_argument("--index", action="store", dest="index", help="Aligner index")
parser.add_argument("--aligner", action="store", dest="aligner", default="bwa",
                    help="Aligner to use. Allowed: bwa(default), bowtie2")
parser.add_argument("-x", "--max_insert_size", action="store", dest="max_insert_size", type=int,
                    help="Maximum insert size to treat read pair as concordant. Affects bowtie2 only.")

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
                     aligner=args.aligner,
                     samples_to_handle=args.samples, threads=args.threads, trimmomatic_dir=args.trimmomatic_dir,
                     trimmer_dir=args.trimmer_dir, bam_util_dir=args.bamutil_dir,
                     mismatch_number=args.mismatch_number, pe_reads_score=args.pe_score,
                     se_read_score=args.se_score, min_adapter_len=args.min_adapter_len,
                     sliding_window_size=args.sliding_window_size,
                     average_quality_threshold=args.average_quality_threshold,
                     leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
                     crop_length=None, head_crop_length=None, min_len=args.min_len,
                     base_quality=args.base_quality,
                     remove_intermediate_files=not args.keep_intermediate_files,
                     filtered_reads=args.filtered,
                     aligned_reads=args.aligned, aligned_and_clipped_reads=args.aligned_and_clipped,
                     max_insert_size=args.max_insert_size,
                     min_coverage_for_filtering=args.min_coverage_for_filtering,
                     max_coverage_for_variant_call=args.max_coverage_for_variant_call,
                     adjust_mapping_quality=args.adjust_mapping_quality)


