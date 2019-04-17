#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_cluster_dir", action="store", dest="input_cluster_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with files with clusters")
parser.add_argument("-s", "--input_seq_dir", action="store", dest="input_seq_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with files with sequences")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory to write_output")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="families",
                    help="Extraction mode. Allowed - 'families', 'species'. Default - 'families'")
parser.add_argument("-w", "--white_list_ids", action="store", dest="white_list_ids",
                    help="File with ids from white list. ")

parser.add_argument("-e", "--seq_extension", action="store", dest="seq_extension", default="fasta",
                    help="Extension of files in sequence directory. Default - 'fasta'")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of files in sequence directory. Default - 'fasta'")
parser.add_argument("-l", "--label_ids", action="store_true", dest="label_ids",
                    help="Label ids by species names. Default - don't label")

parser.add_argument("-g", "--separator_for_labeling", action="store", dest="separator_for_labeling", default="@",
                    help="Separator to use for labeling. Default - '@'")
parser.add_argument("-r", "--label_last", action="store_false", dest="label_first", default=True,
                    help="Place label at the end of id")

args = parser.parse_args()

FileRoutines.safe_mkdir(args.output_dir)
SequenceClusterRoutines.extract_sequences_by_clusters(args.input_cluster_dir, args.input_seq_dir, args.output_dir,
                                                      file_with_white_list_cluster_ids=args.white_list_ids,
                                                      mode=args.mode, sequence_file_extension=args.seq_extension,
                                                      sequence_file_format=args.format,
                                                      label_species=args.label_ids,
                                                      separator_for_labeling=args.separator_for_labeling,
                                                      species_label_first=args.label_first)
