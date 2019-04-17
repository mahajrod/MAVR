#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import shutil
from Bio import SeqIO
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file.")
parser.add_argument("-n", "--number_of_sequences", action="store", dest="number_of_sequences", default=5000, type=int,
                    help="Number of sequences per splited file. Default: 5000")
parser.add_argument("-s", "--splited_directory", action="store", dest="splited_directory", default="splited_dir",
                    help="Directory with splited files. Default: 'splited_dir' ")
parser.add_argument("-d", "--output_directory", action="store", dest="output_directory", default="output_dir",
                    help="Directory with output files. Default: 'output_dir'")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Output file prefix")
parser.add_argument("-c", "--command", action="store", dest="command", required=True,
                    help="Command to parallel with xargs")
parser.add_argument("-m", "--common_options", action="store", dest="common_options",
                    help="Options common to all threads")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads")
parser.add_argument("-r", "--retain_temp_dirs", action="store_true", dest="retain_temp_dirs", default=False,
                    help="Store temporary directories")
args = parser.parse_args()

sequence_dict = SeqIO.index_db("temp.idx", args.input, format="fasta")

try:
    os.mkdir(args.splited_directory)
except OSError:
    pass

try:
    os.mkdir(args.output_directory)
except OSError:
    pass

split_index = 1
records_written = 0
record_ids_list = list(sequence_dict.keys())
number_of_records = len(record_ids_list)


while (records_written + args.number_of_sequences) <= number_of_records:

    SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict,
                                                        record_ids_list[records_written:records_written+args.number_of_sequences],
                                                        verbose=True),
                "%s/%s_%i.fasta" % (args.splited_directory, args.output_prefix, split_index), format="fasta")
    split_index += 1
    records_written += args.number_of_sequences

if records_written != number_of_records:
    SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict,
                                                        record_ids_list[records_written:],
                                                        verbose=True),
                "%s/%s_%i.fasta" % (args.splited_directory, args.output_prefix, split_index), format="fasta")

os.remove("temp.idx")

xargs_string = "ls %s | xargs -P %i -I NAME %s -query %s/NAME -out %s/NAME.hits %s" % (args.splited_directory, args.threads,
                                                                                       args.command, args.splited_directory,
                                                                                       args.output_directory, args.common_options)
os.system(xargs_string)

merge_string = "cat %s/* > %s.hits" % (args.output_directory, args.output_prefix)

os.system(merge_string)

if not args.retain_temp_dirs:
    shutil.rmtree(args.output_directory)
    shutil.rmtree(args.splited_directory)

