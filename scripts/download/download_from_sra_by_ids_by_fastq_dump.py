#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Tools.NCBIToolkit import FastqDump


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--ids", action="store", dest="ids",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of SRA ids to download")
parser.add_argument("-f", "--id_file", action="store", dest="id_file",
                    help="File with SRA ids(one per line) to download")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of simultaneous downloads")
parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", default="./",
                    help="Output directory. Default: current directory")

args = parser.parse_args()

if (not args.ids) and (not args.id_file):
    raise ValueError("Both ids and id file were not set")

id_list = IdList(filename=args.id_file) if args.id_file else args.ids

FastqDump.threads = args.threads
FastqDump.parallel_download(id_list, args.out_dir, split_pe=True, retain_original_ids=True)
