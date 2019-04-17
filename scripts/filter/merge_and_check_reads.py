#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from multiprocessing import Pool
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="sample_dir",
                    help="Directory with samples")
parser.add_argument("-a", "--samples", action="store", dest="samples",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples."
                         "In sample directory should one(in case SE reads) or two(in case PE reads) files."
                         "Filenames should should contain '_1.fq' or '_1.fastq' for forward(left) reads, "
                         " '_2.fq' or '_2.fastq' for reverse(right) reads and '.fq' or '.fastq' for SE reads."
                         "Also files could be compressed by gzip")

parser.add_argument("-z", "--gziped_input", action="store_true", dest="gziped_input",
                    help="Files with reads are compressed by gzip")
parser.add_argument("-c", "--merging_threads", action="store", dest="merging_threads", type=int,
                    default=1,
                    help="Number of threads to use during merging")

parser.add_argument("-r", "--merged_dir", action="store", dest="merged_dir",
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    default="./fastq/", help="Directory to write merged files. Default: ./fastq/")
parser.add_argument("-k", "--kmer_dir", action="store", dest="merged_dir",
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    default="./kmers/", help="Directory to write results of kmer analysis. Default: ./kmers/")

parser.add_argument("-i", "--input_file", action="store", dest="input", type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of fasta or fastq files or directories containing them.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["svg", "png", "jpg"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default -10")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=int, default=1000000,
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L)")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use by jellyfish. Default - 1")
parser.add_argument("-b", "--count_both_strands", action="store_true", dest="count_both_strands",
                    help="Count kmers in both strands. NOTICE: only mer or its reverse-complement, whichever "
                         "comes first lexicographically, is stored and the count value is the number of "
                         "occurrences of both. So this option is not suitable for generating sets of forward "
                         "and reverse-complement kmers. For this case use -r/--add_reverse_complement option. "
                         "Not compatible with -r/--add_reverse_complement option.")
parser.add_argument("-r", "--add_reverse_complement", action="store_true", dest="add_rev_com",
                    help="Add reverse-complement sequences before counting kmers. "
                         "Works only for fasta sequences. "
                         "Not compatible with -b/--count_both_strands option")
parser.add_argument("-j", "--jellyfish_path", action="store", dest="jellyfish_path",
                    help="Path to jellyfish")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=5,
                    help="Low limit of histogram without logscale")
parser.add_argument("-g", "-high_limit", action="store", dest="high_limit", type=int, default=100,
                    help="High limit of histogram without logscale")


args = parser.parse_args()

samples = args.samples if args.samples else os.listdir(args.samples_dir)

merging_pool = Pool(args.merging_threads)

# TODO: write or delete this script