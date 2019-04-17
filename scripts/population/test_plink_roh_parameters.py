#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Population import PLINK

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf files")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output direct with sequences")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-d", "--plink_bin_dir", action="store", dest="plink_bin_dir", default="",
                    help="Directory with PLINK binary")


args = parser.parse_args()

PLINK.path = args.plink_bin_dir
PLINK.test_roh_parameters(args.output_dir,
                          args.output_prefix,
                          args.input_vcf,
                          allow_noncanonical_chromosome_names=True,
                          keep_autoconverted_files=None,
                          window_length_in_kb=None,
                          min_homozygous_snps_per_window=(2, 51, 1),
                          min_homozygous_snps_in_roh=(2, 101, 1),
                          max_heterozygous_snps_per_window=(1, 11, 1),
                          max_heterozygous_snps=(1, 21, 1),
                          max_inverse_density_of_homozygous_snps_in_kb_per_snp=(50, 1000, 50),
                          )
