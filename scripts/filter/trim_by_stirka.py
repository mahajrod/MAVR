#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Stirka import Trimmer


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True, type=Trimmer.check_dir_path,
                    help="Input directory with samples. Files of each sample have to be located in its own dir")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                    help="Comma-separated list of samples to trim. Reads in sample directories have to be named in "
                         "following manner: <Sample_name>_1.fastq, <Sample_name>_2.fastq. If not set all directories in"
                         "input directory will be treated as sample dir")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True, type=Trimmer.check_dir_path,
                    help="Input directory with samples. Files of each sample have to be located in its own dir")
parser.add_argument("-a", "--adapters", action="store", dest="adapters", required=True,
                    help="File with adapter kmers")

parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use. Default: 1")
parser.add_argument("-d", "--stirka_bin_dir", action="store", dest="stirka_bin_dir", default="",
                    help="Directory with trimmer binary")
parser.add_argument("-w", "--write_fasta", action="store_true", dest="write_fasta", default=False,
                    help="Write fasta file. Default: False")
parser.add_argument("-f", "--format", action="store", dest="format", default="fastq",
                    help="Format of input files. Allowed: reads, fasta, fastq(default)")

args = parser.parse_args()

Trimmer.threads = args.threads
Trimmer.path = args.stirka_bin_dir
Trimmer.mass_trim(args.input_dir, args.output_dir, args.adapters, sample_list=args.sample_list,
                  write_fasta=args.write_fasta, input_format=args.format)
