#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import pandas as pd
import sys
import argparse
from RouToolPa.Routines import FilteringRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-k", "--kraken_output", action="store", dest="kraken_output", required=True,
                    help="File with kraken output. Might be gzipped")
parser.add_argument("-t", "--taxon_id_list", action="store", dest="taxon_id_list", required=True,
                    type=lambda s: pd.read_csv(s, header=None, dtype=str).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of taxon ids")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file to write extracted sequences. Default: stdout")


args = parser.parse_args()
print(args.taxon_id_list)
print(type(args.taxon_id_list[0]))
FilteringRoutines.extract_seq_ids_by_taxa_from_kraken_report(args.kraken_output, args.taxon_id_list, args.output)
