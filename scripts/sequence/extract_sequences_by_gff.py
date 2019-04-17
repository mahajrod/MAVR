#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import AnnotationsRoutines
#from RouToolPa.Routines.Sequence import record_generator


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-t", "--type", action="store", dest="type", default="gene", type=lambda s: s.split(","),
                    help="Comma-separated list of feature types to extract")
parser.add_argument("-g", "--gff_file", action="store", dest="gff_file",
                    help="Gff file with annotations to extract")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db', 'index'(default), 'parse'")

args = parser.parse_args()


AnnotationsRoutines.extract_sequences_by_gff(args.input, args.gff_file, args.output, type_list=args.type,
                                             parsing_mode=args.parsing_mode, tmp_index_file="temp.idx",
                                             format=args.format)


"""
tmp_index_file = "temp.idx"
args.type = args.type.split(",")


annotations_dict = SeqIO.to_dict(GFF.parse(open(args.gff_file)))
print annotations_dict
print("Parsing %s..." % args.input)
sequence_dict = SequenceRoutines.parse_seq_file(args.input, args.parsing_mode, args.format, index_file=tmp_index_file ) # SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)


SeqIO.write(SequenceRoutines.record_generator(annotations_dict, sequence_dict, args.type), args.output, format=args.format)

if args.parsing_mode == "index_db":
    os.remove(tmp_index_file)
"""



