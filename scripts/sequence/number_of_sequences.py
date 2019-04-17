#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from Bio import SeqIO
from RouToolPa.Routines.File import make_list_of_path_to_files, detect_filetype_by_extension, split_filename



def file_filter(filename):
    allowed_filetypes = ["fasta", "genbank"]
    return True if detect_filetype_by_extension(filename) in allowed_filetypes else False

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output",  default="stdout",
                    help="Output file with number of sequences in file(files). Default: stdout")
parser.add_argument("-a", "--write_header", action="store_true", dest="write_header",
                    help="Write header in output file. Default: false")
parser.add_argument("-i", "--input", action="store", dest="input", type=lambda s: s.split(","),
                    help="Comma-separated list of files with sequences or directories containing them",  required=True)
parser.add_argument("-w", "--write_dir_path", action="store_true", dest="write_dir_path",
                    help="Write directory name(if directory is source of files) in output file. Default: false")
parser.add_argument("-e", "--write_ext", action="store_true", dest="write_ext",
                    help="Write extensions of files with sequences in output file. Default: false")
args = parser.parse_args()

files_list = sorted(make_list_of_path_to_files(args.input, file_filter))


out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

if args.write_header:
    out_fd.write("#file/sample\tnumber_of_sequences\n")
for filename in files_list:
    if args.output != "stdout":
        print("Counting variants in %s ..." % filename)
    directory, prefix, extension = split_filename(filename)
    filetype = detect_filetype_by_extension(filename)
    number_of_sequences = 0
    with open(filename, "r") as seq_fd:
        try:
            for record in SeqIO.parse(seq_fd, filetype):
                number_of_sequences += 1
        except:
            if args.output != "stdout":
                print("%s was not opened" % filename)

    if args.write_dir_path and args.write_ext:
        name = filename
    elif args.write_dir_path:
        name = (directory + prefix) if directory else prefix
    elif args.write_ext:
        name = prefix + extension
    else:
        name = prefix

    out_fd.write("%s\t%i\n" % (name, number_of_sequences))

if args.output != "stdout":
    out_fd.close()