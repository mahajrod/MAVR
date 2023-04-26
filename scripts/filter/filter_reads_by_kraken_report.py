#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import pd
import argparse
from RouToolPa.Routines import FilteringRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_fastq", action="store", dest="forward_fastq", required=True,
                    help="Forward fastq file. Might be gzipped.")
parser.add_argument("-r", "--reverse_fastq", action="store", dest="reverse_fastq",
                    help="Reverse fastq file. Might be gzipped.")
parser.add_argument("-k", "--kraken_output", action="store", dest="kraken_output", required=True,
                    help="File with kraken output. Might be gzipped")
parser.add_argument("-t", "--taxon_id_list", action="store", dest="taxon_id_list", required=True,
                    type=lambda s: list(pd.read_csv(s, header=None, dtype=str).squeeze("columns")) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of taxon ids")
parser.add_argument("-c", "--check_read_id", action="store_true", dest="check_read_id", default=False,
                    help="Check correspondence of read ids in read files and KRAKEN output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")


args = parser.parse_args()

FilteringRoutines.extract_reads_by_kraken_report(args.forward_fastq, args.kraken_output, args.taxon_id_list,
                                                 args.output_prefix, reverse_fastq=args.reverse_fastq, gzip=True,
                                                 check_read_id=args.check_read_id)
