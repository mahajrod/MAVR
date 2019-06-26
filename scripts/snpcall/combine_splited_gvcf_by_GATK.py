#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import CatVariants


parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gvcf file")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-i", "--input_gvcf_list", action="store", dest="input_gvcf_list",
                    type=CatVariants.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of gvcf files to combine",  required=True,)
parser.add_argument("-s", "--sorted", action="store_true", dest="sorted", default=False,
                    help="Input gvcf are coordinate sorted. Default: False")
parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-m", "--remove_intermediate_files", action="store_true", dest="remove_intermediate_files", default=False,
                    help="Remove intermediate files. Default: False")
parser.add_argument("-d", "--tmp_dir", action="store", dest="tmp_dir", default="./tmp_combine_gvcf/",
                    help="Directory for temporary files. Default: ./tmp_combine_gvcf/")
parser.add_argument("-x", "--max_files_per_merging", action="store", dest="max_files_per_merging", default=50, type=int,
                    help="Maximum number of files per merging. Default: 50")

args = parser.parse_args()

CatVariants.jar_path = args.gatk_dir
CatVariants.threads = args.threads
#print args.extension_list
CatVariants.combine_gvcf(args.reference, args.input_gvcf_list, args.output,
                         input_is_sorted=args.sorted,
                         extension_list=args.extension_list,
                         tmp_dir=args.tmp_dir,
                         max_files_per_merging=args.max_files_per_merging,
                         iteration=0,
                         threads=None,
                         remove_intermediate_files=args.remove_intermediate_files)
