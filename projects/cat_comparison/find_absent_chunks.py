__author__ = 'mahajrod'
import os
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_file", required=True,
                    help="Input dir with exonerate output named in style *_<chunk_number>.out")
parser.add_argument("-t", "--total_number_of_chunks", action="store", dest="number_of_chunks", type=int,
                    required=True,
                    help="Total number of chunks")

args = parser.parse_args()


chunk_files = os.listdir(args.input_dir)

chunk_numbers = sorted([int(n.split(".")[-2].split("_")[-1]) for n in chunk_files])

if chunk_numbers[-1] > args.number_of_chunks:
    print("Largest number of present chunks(%i) is larger than expected(%i)" % (chunk_numbers[-1],
                                                                                args.number_of_chunks))
print("Absent chunks:")

for i in range(1, args.number_of_chunks + 1):
    if i not in chunk_numbers:
        print i
