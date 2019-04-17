#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import argparse
from RouToolPa.Routines import FileRoutines
from RouToolPa.Collections.General import SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_files_dir", action="store", dest="input_files_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with files names by taxa id")
parser.add_argument("-o", "--output_files_dir", action="store", dest="output_files_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory to write fam files named by species names")
parser.add_argument("-d", "--syn_file", action="store", dest="syn_file", required=True,
                    help="File with taxa ids and species names")
parser.add_argument("-k", "--key_index", action="store", dest="key_index", type=int, default=0,
                    help="Key column in file with synonyms(0-based). Default: 0")
parser.add_argument("-v", "--value_index", action="store", dest="value_index", type=int, default=1,
                    help="Value column in file with synonyms(0-based). Default: 1")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Prefix of comments in synonyms file. Default - '#'")
parser.add_argument("-m", "--columns_separator", action="store", dest="separator", default="\t",
                    help="Column separator in file with synonyms")
parser.add_argument("-e", "--header", action="store_true", dest="header", default=False,
                    help="Header is present in synonyms file. Default - False")

args = parser.parse_args()

syn_dict = SynDict()
syn_dict.read(args.syn_file, header=args.header, separator=args.separator, key_index=args.key_index,
              value_index=args.value_index, comments_prefix=args.comments_prefix)

FileRoutines.safe_mkdir(args.output_files_dir)
input_files = os.listdir(args.input_files_dir)
for filename in input_files:
    directory, taxon_id, extension = FileRoutines.split_filename(filename)
    if taxon_id not in syn_dict:
        print("Species name was not found for taxon %s" % taxon_id)
        continue
    shutil.copy("%s%s" % (args.input_files_dir, filename),
                "%s%s%s" % (args.output_files_dir, syn_dict[taxon_id], extension))
