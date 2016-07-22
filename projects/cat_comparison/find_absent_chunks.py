#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input dir with exonerate output named in style *_<chunk_number>.out")
parser.add_argument("-n", "--total_number_of_chunks", action="store", dest="number_of_chunks", type=int,
                    required=True,
                    help="Total number of chunks")

args = parser.parse_args()


chunk_files = os.listdir(args.input)

chunk_numbers = sorted([int(n.split(".")[-2].split("_")[-1]) for n in chunk_files])

if chunk_numbers[-1] > args.number_of_chunks:
    print("Largest number of present chunks(%i) is larger than expected(%i)" % (chunk_numbers[-1],
                                                                                args.number_of_chunks))

absent_chunks = []
for i in range(1, args.number_of_chunks + 1):
    if i not in chunk_numbers:
        absent_chunks.append(i)

if absent_chunks:
    print("Absent chunks(total %i): " % len(absent_chunks))
    for k in absent_chunks:
        print(k)
else:
    print("No absent chunks")
