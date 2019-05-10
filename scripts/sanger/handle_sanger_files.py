#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from Pipelines import SangerPipeline


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with Sanger (.ab1) files")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", default="./",
                    help="LAST database")
parser.add_argument("-p", "--o utput_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--read_subfolders", action="store_true", dest="read_subfolders", default=False,
                    help="Read data from subfolders of input directory")
parser.add_argument("-l", "--min_length", action="store", dest="min_length", default=50, type=int,
                    help="Minimum length of trimmed sequence to retain. Default: 50")
parser.add_argument("-n", "--min_median_qual", action="store", dest="min_median_qual", default=15, type=int,
                    help="Minimum median quality of sequence to retain. "
                         "Can be applied along with -k/--min_mean_qual option. Default: 15")
parser.add_argument("-k", "--min_mean_qual", action="store", dest="min_mean_qual", default=15, type=int,
                    help="Minimum mean quality of sequence to retain. "
                         "Can be applied along with -n/--min_median_qual option. Default: 15")


args = parser.parse_args()

SangerPipeline.workdir = args.output_dir
SangerPipeline.handle_sanger_data(args.input_dir, args.output_prefix, outdir=None,
                                  read_subfolders=args.read_subfolders,
                                  min_mean_qual=args.min_mean_qual, min_median_qual=args.min_median_qual,
                                  min_len=args.min_length)
