#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.SRAToolkit import FastqDump


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with sra files each sra file must be located in directory named by id")

parser.add_argument("-s", "--sra_id_list", action="store", dest="sra_id_list", type=lambda s: s.split(","),
                    help="Comma-separated list of sra ids to extract. If not set all folder_names in "
                         "input directory will be treated as ids")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads")
parser.add_argument("-u", "--unpaired", action="store_false", dest="paired", default=True,
                    help="Sra archives contain unpaired reads")
parser.add_argument("-p", "--fastq_dump_dir", action="store", dest="fastq_dump_dir", default="",
                    help="Path to directory with fastq-dump directory")
parser.add_argument("-d", "--dont_retain_original_ids", action="store_false", dest="retain_original_ids", default=True,
                    help="Don't retain original read ids, use sra ids instead")

args = parser.parse_args()

FastqDump.threads = args.threads
FastqDump.path = args.fastq_dump_dir
FastqDump.parallel_unpack(args.input_dir, args.output_dir, sra_id_list=args.sra_id_list, paired_input=args.paired,
                          retain_original_ids=args.retain_original_ids)
