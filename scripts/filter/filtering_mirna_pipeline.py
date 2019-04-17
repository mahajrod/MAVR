#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Pipelines import FilteringPipeline
from RouToolPa.Routines.File import check_path


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")

#parser.add_argument("-x", "--general_stat_file", action="store", dest="general_stat_file", required=True,
#                    help="File to write general statistics about filtration")

parser.add_argument("-s", "--samples", action="store", dest="samples", type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-a", "--adapters", action="store", dest="adapters",
                    required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of adapters or files with adapters to trim by Cutadapt")
parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=20,
                    type=int,
                    help="Quality threshold for sliding window or whole read."
                         "Depends on -q/--average_quality_threshold option.Default - 20.")
parser.add_argument("-b", "--base_quality", action="store", dest="base_quality", default="phred33",
                    help="Type of base quality. Possible variants: phred33, phred64. Default - phred33 ")
parser.add_argument("-l", "--min_length", action="store", dest="min_len", type=int, default=15,
                    help="Minimum length of read to retain. Default - 15")
parser.add_argument("-e", "--max_length", action="store", dest="max_len", type=int, default=39,
                    help="Maximum length of read to retain. Default - 39")

parser.add_argument("-c", "--cutadapt_dir", action="store", dest="cutadapt_dir", default="",
                    help="Path to Cutadapt directory")
parser.add_argument("-f", "--facut_dir", action="store", dest="facut_dir", default="",
                    help="Path to Facut directory")
parser.add_argument("-r", "--remove_intermediate_files", action="store_true",
                    dest="remove_intermediate_files", default=False,
                    help="Remove intermediate files")
parser.add_argument("-z", "--read_name_type", action="store", dest="read_name_type", default="illumina",
                    help="Read name type")


args = parser.parse_args()

"""
EXAMPLE
skliver@supermicro:
cd ~/workdir/yeast/nizhnikov/good_run/fastq
~/soft/MAVR/scripts/filter/filtering_pipeline.py -d raw/ -o ./ \
                                                 -k ~/data/service_seq/trueseq_adapters_with_rev_com_23_mer.kmer
                                                 -a ~/soft/Trimmomatic-0.35/adapters/TruSeq3-PE.fa
                                                 -j ~/soft/Trimmomatic-0.35/
                                                 -x filtering_general.stat -r
"""


FilteringPipeline.filter_mirna(args.samples_dir, args.output_dir, args.adapters,
                               general_stat_file=None,
                               samples_to_handle=args.samples, cutadapt_dir=args.cutadapt_dir,
                               min_len=args.min_len, max_len=args.max_len,
                               remove_intermediate_files=args.remove_intermediate_files,
                               average_quality_threshold=args.average_quality_threshold, base_quality=args.base_quality, read_name_type=args.read_name_type
                               )

