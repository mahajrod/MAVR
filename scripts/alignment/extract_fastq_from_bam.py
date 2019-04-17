#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Samtools import SamtoolsV1
from RouToolPa.Tools.Bedtools import BamToFastq
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input bam file")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use. Default - 1")
parser.add_argument("-p", "--prepare_bam", action="store_true", dest="prepare_bam",
                    help="Prepare bam for reads extraction(filter out supplementary and not primary alignments"
                         "and sort by name)")
"""
parser.add_argument("-e", "--prepared_bam", action="store", dest="prepared_bam",
                    help="File to write sorted bam file. Required if -p/--prepare_bam option is set")
"""
parser.add_argument("-e", "--prepared_bam_prefix", action="store", dest="prepared_bam_prefix",
                    help="Prefix of sorted bam file(s). Required if -p/--prepare_bam option is set")
parser.add_argument("-d", "--temp_dir", action="store", dest="temp_dir",
                    help="Directory to use for temporary files. Required if -p/--prepare_bam option is set")
parser.add_argument("-o", "--out_prefix", action="store", dest="out_prefix", required=True,
                    help="Prefix of output fastq files")
parser.add_argument("-s", "--single_ends", action="store_false", dest="paired", default=True,
                    help="Reads are SE")
parser.add_argument("-x", "--mix_ends", action="store_true", dest="mix_ends", default=False,
                    help="Reads are mix of PE and SE")
parser.add_argument("-m", "--max_memory_per_thread", action="store", dest="max_memory_per_thread", default="1G",
                    help="Maximum memory per thread. Default - 1G")
args = parser.parse_args()

if args.prepare_bam and ((not args.prepared_bam_prefix) or (not args.temp_dir)):
    raise ValueError("Options -e/--prepared_bam_prefix and -m/--temp_dir must be set if -p/--prepare_bam option is used")

SamtoolsV1.threads = args.threads

if args.prepare_bam or args.mix_ends:
    FileRoutines.safe_mkdir(FileRoutines.check_path(args.temp_dir))
    prepared_pe_bam_file = "%s.bam" % args.prepared_bam_prefix
    prepared_unpaired_bam_file = ("%s.unpaired.bam" % args.prepared_bam_prefix) if args.mix_ends else None
    """
    SamtoolsV1.prepare_bam_for_read_extraction(args.input, args.prepared_bam, temp_file_prefix=args.temp_dir,
                                               max_memory_per_thread=args.max_memory_per_thread)
    """
    SamtoolsV1.prepare_bam_for_read_extraction(args.input, prepared_pe_bam_file, temp_file_prefix=args.temp_dir,
                                               max_memory_per_thread=args.max_memory_per_thread,
                                               bam_file_to_write_unpaired_reads=prepared_unpaired_bam_file)
if args.paired:
    left_fastq = "%s_1.fastq" % args.out_prefix
    right_fastq = "%s_2.fastq" % args.out_prefix
    unpaired_fastq = "%s.unpaired.fastq" % args.out_prefix
else:
    left_fastq = "%s.fastq" % args.out_prefix
    right_fastq = None

if args.mix_ends:
    BamToFastq.convert(prepared_unpaired_bam_file, unpaired_fastq, out_right_fastq=None)

#BamToFastq.convert(args.prepared_bam if args.prepare_bam else args.input, left_fastq, out_right_fastq=right_fastq)
BamToFastq.convert(prepared_pe_bam_file if args.prepare_bam else args.input, left_fastq, out_right_fastq=right_fastq)
