#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import GenotypeGVCFs

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--splited_prefix", action="store", dest="splited_prefix", required=True,
                    help="Prefix of splited output files")
parser.add_argument("-d", "--splited_dir", action="store", dest="splited_dir", required=True,
                    help="Directory to write splited vcf files")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output vcf files with all regions")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-p", "--picard_directory", action="store", dest="picard_dir", default="",
                    help="Directory with Picard jar")
parser.add_argument("-i", "--gvcf_list", action="store", dest="gvcf_list",
                    type=GenotypeGVCFs.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of gvcf files to genotype",  required=True,)
parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-m", "--max_alternate_alleles", action="store", dest="max_alternate_alleles", type=int,
                    help="Maximum number of alternative allels. Default: GATK default")

args = parser.parse_args()

GenotypeGVCFs.jar_path = args.gatk_dir
GenotypeGVCFs.threads = args.threads
GenotypeGVCFs.parallel_genotype(args.reference, args.gvcf_list, args.splited_dir, args.splited_prefix, args.output,
                                max_total_scaffold_length_per_chunk=100000,
                                max_scaffold_number_per_chunk=5, length_dict=None,
                                parsing_mode="parse", region_list=None,
                                extension_list=args.extension_list,
                                max_alternate_alleles=args.max_alternate_alleles,
                                picard_jar_path=args.picard_dir)
