#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from RouToolPa.Routines import TreeFamRoutines, FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--families_id_file", action="store", dest="fam_id_file",
                    help="File with ids of families to extract. Extract all families if not set")
parser.add_argument("-f", "--families_file", action="store", dest="fam_file", required=True,
                    help="File with families")
parser.add_argument("-p", "--pep_file", action="store", dest="pep_file", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="List of comma-separated files/directories with proteins")
parser.add_argument("-r", "--pep_file_format", action="store", dest="pep_file_format", default="fasta",
                    help="Format of file with proteins")
parser.add_argument("-c", "--create_dir_for_each_family", action="store_true", dest="create_dir_for_each_family",
                    help="Create separate directory for each family")
parser.add_argument("-o", "--output_prefix", action="store", dest="output",
                    help="Output prefix to use")
parser.add_argument("-d", "--output_directory", action="store", dest="out_dir", default="./",
                    help="Directory to write output")

args = parser.parse_args()

TreeFamRoutines.extract_proteins_from_selected_families(args.fam_id_file, args.fam_file, args.pep_file,
                                                        output_dir=args.out_dir, pep_format=args.pep_file_format,
                                                        out_prefix=args.output,
                                                        create_dir_for_each_family=args.create_dir_for_each_family)
