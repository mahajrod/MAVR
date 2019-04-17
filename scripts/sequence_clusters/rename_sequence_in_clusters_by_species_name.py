#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Routines import SequenceRoutines, FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with per family sequences labeled by  species")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory to write output")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator between species name and sequence id. Default - '@'")
parser.add_argument("-r", "--label_last", action="store_false", dest="label_first", default=True,
                    help="Species label is at the end of id")

args = parser.parse_args()

input_file_list = sorted(os.listdir(args.input_dir))
FileRoutines.safe_mkdir(args.output_dir)
if args.label_first:
    id_expression = lambda record_id: record_id.split(args.separator)[0]
else:
    id_expression = lambda record_id: record_id.split(args.separator)[1]


for filename in input_file_list:
    input_file = "%s%s" % (args.input_dir, filename)
    output_file = "%s%s" % (args.output_dir, filename)
    SequenceRoutines.rename_records_from_files(input_file, output_file, record_id_expression=id_expression,
                                               clear_description=True)



