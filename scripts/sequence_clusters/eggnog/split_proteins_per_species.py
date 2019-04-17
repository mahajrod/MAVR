#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import EggNOGRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_pep_dir", action="store", dest="input_pep_dir", required=True,
                    type=EggNOGRoutines.check_path,
                    help="Directory with input proteins(splited by families or something like)")
parser.add_argument("-o", "--output_pep_dir", action="store", dest="output_pep_dir", required=True,
                    type=EggNOGRoutines.check_path,
                    help="Directory to write proteins splited by species")

args = parser.parse_args()

EggNOGRoutines.split_proteins_per_species(args.input_pep_dir, args.output_pep_dir,
                                          input_format="fasta", output_format="fasta")
