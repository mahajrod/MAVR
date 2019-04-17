#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Emapper


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--eggnog_db", action="store", dest="eggnog_db", required=True,
                    help="EggNOG database to use for ortholog assignment")
parser.add_argument("-e", "--eggnog_db_prefix", action="store", dest="eggnog_db_prefix",
                    help="Prefix of EggNOG database ids")
parser.add_argument("-s", "--species", action="store", dest="species", required=True,
                    help="Species name. Required for labels in fam_file")
parser.add_argument("-l", "--label_separator", action="store", dest="label_separator", default="@",
                    help="Label separator.Default: @")
parser.add_argument("-i", "--input_seq", action="store", dest="input",
                    help="Input file with sequences")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("--emapper_dir", action="store", dest="path", default="",
                    help="Path to directory with emapper.py")
parser.add_argument("--temp_dir", action="store", dest="temp_dir",
                    help="Path to temp directory")

args = parser.parse_args()

"""
Usage:
assign_orthologs_emapper.py -d veNOG -o amazona_vittata.longest_pep -t 30 -i amazona_vittata.longest_pep.pep
"""

Emapper.threads = args.threads
Emapper.path = args.path
Emapper.timelog = None
Emapper.tmp_dir = args.temp_dir

Emapper.assign_orthologs(args.input, args.eggnog_db, args.output_prefix,
                         db_type=None, query_type=None, target_orthologs=None,
                         resume=None, report_orthologs=True, override_last_run=None,
                         no_comments_in_output=None, no_refine=None, no_annotation=None,
                         no_search=None, usemem=True, eggnogdb_prefix=args.eggnog_db_prefix,
                         species_name=args.species, label_separator=args.label_separator)

