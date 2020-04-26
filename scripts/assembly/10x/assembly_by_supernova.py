#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import matplotlib

from RouToolPa.Tools.Assemblers import Supernova

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--fastq_dir", action="store", dest="fastq_dir", required=True,
                    help="Directory with fastqs")
parser.add_argument("-t", "--threads", action="store", dest="threads", required=True, type=int,
                    help="Number of threads to use")
parser.add_argument("-m", "--max_memory", action="store", dest="max_memory", required=True, type=int,
                    help="Maximum memory to use (in Gigabytes)")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory(only numbers, letters, dash, and underscore allowed)")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")

args = parser.parse_args()

Supernova.max_memory = args.max_memory
Supernova.threads = args.threads
Supernova.assembly(args.fastq_dir, args.output_dir, args.output_prefix, max_reads="all", disable_ui=True)
