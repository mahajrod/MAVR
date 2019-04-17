#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines import SequenceRoutines, FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of  genbank files/directories with transcript annotations")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with id of transcripts to deal with")

args = parser.parse_args()

if args.id_file:
    id_list = IdList()
    id_list.read(args.id_file)
else:
    id_list = None
SequenceRoutines.extract_introns_from_transcripts_from_genbank_files(args.input, args.output,
                                                                     transcript_id_white_list=id_list)
