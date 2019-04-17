#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences.")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with sequences.")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-n", "--min_length", action="store", dest="min_length", type=int,
                    help="Minimun length of sequence to store. Default: filter not used")
parser.add_argument("-x", "--max_length", action="store", dest="max_length", type=int,
                    help="Maximum length of sequence to store. Default: filter not set")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File to write ids of extracted sequences. Default - don't write")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default='parse',
                    help="Parsing mode. Allowed: parse(default), index, index_db")

args = parser.parse_args()


SequenceRoutines.extract_sequences_by_length_from_file(args.input_file, args.output_file,
                                                       min_len=args.min_length,
                                                       max_len=args.max_length, format=args.format,
                                                       tmp_index_file="tmp.idx",
                                                       id_file=args.id_file,
                                                       parsing_mode=args.parsing_mode)

"""
if (args.min_length is None) and (args.max_length is None):
    raise ValueError("Both minimum and maximum lengths were not set")
elif (args.min_length is not None) and (args.max_length is not None) and (args.min_length > args.max_length):
    raise ValueError("Minimum length is greater then maximum lengths")

tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)

if (args.min_length is not None) and (args.max_length is not None):
    length_expression = lambda record: args.min_length <= len(record.seq) <= args.max_length
elif args.min_length is not None:
    length_expression = lambda record: len(record.seq) >= args.min_length
else:
    length_expression = lambda record: len(record.seq) <= args.max_length

#id_list = filter_sequences(sequence_dict, expression)

SeqIO.write(SequenceRoutines.record_by_expression_generator(sequence_dict, expression=length_expression,
                                                            id_file=args.id_file),
            args.output_file, format=args.format)

#SeqIO.write(record_by_id_generator(sequence_dict, id_list), args.output_file, format=args.format)

os.remove(tmp_index_file)
"""






