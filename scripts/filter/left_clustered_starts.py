#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from collections import OrderedDict
from RouToolPa.Routines.Functions import get_cigar_str_len



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")
parser.add_argument("-s", "--min_cluster_size", action="store", dest="min_cluster_size", default=5, type=int,
                    help="Minimum size of clustered starts")
parser.add_argument("-f", "--histo_file", action="store", dest="histo_file", default="histo.t",
                    help="File with histogram")

args = parser.parse_args()

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

read_counter = 0


for line in in_fd:
    if line[0] != "@":
        tmp = line
        break
    else:
        out_fd.write(line)


tmp_list = tmp.strip().split()
flag = int(tmp_list[1])
cigar_string = tmp_list[5]
length = get_cigar_str_len(cigar_string)
strand = -1 if flag & 16 == 16 else 1
position = int(tmp_list[3])
alignment_start = position if strand == 1 else position + length - 1

start_dict = {1: OrderedDict({}), -1: OrderedDict({})}
histo_dict = {1: OrderedDict({}), -1: OrderedDict({})}

start_dict[strand][alignment_start] = [tmp]
read_counter += 1

for line in in_fd:
    tmp_list = line.strip().split()
    flag = int(tmp_list[1])
    cigar_string = tmp_list[5]
    length = get_cigar_str_len(cigar_string)
    strand = -1 if flag & 16 == 16 else 1
    position = int(tmp_list[3])
    alignment_start = position if strand == 1 else position + length - 1
    if alignment_start in start_dict[strand]:
        start_dict[strand][alignment_start].append(line)

    else:
        start_dict[strand][alignment_start] = [line]

    read_counter += 1

    if read_counter == 10000:
        read_counter = 0
        for strand in 1, -1:
            for key in start_dict[strand]:
                if position <= key:
                    break
                count = len(start_dict[strand][key])
                if count >= args.min_cluster_size:
                    histo_dict[strand][key] = count
                    for string in start_dict[strand][key]:

                        out_fd.write(string)

                del start_dict[strand][key]
#print (start_dict)

for strand in 1, -1:
    for key in start_dict[strand]:
        #print (len(start_dict[strand][key]))
        count = len(start_dict[strand][key])
        if count >= args.min_cluster_size:
            histo_dict[strand][key] = count
            for string in start_dict[strand][key]:
                out_fd.write(string)

        del start_dict[strand][key]

with open(args.histo_file, "w") as histo_fd:
    histo_fd.write("#strand\tstart\tcoverage\n")
    for strand in 1, -1:
        for key in histo_dict[strand]:
            histo_fd.write("%i\t%i\t%i\n" % (strand, key, histo_dict[strand][key]))

if args.output != "output":
    out_fd.close()
if args.input != "stdin":
    in_fd.close()

