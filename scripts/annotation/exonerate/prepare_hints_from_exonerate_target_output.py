#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Exonerate, AUGUSTUS
from RouToolPa.Routines.File import make_list_of_path_to_files


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Files with exonerate target output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-d", "--white_id_file", action="store", dest="white_id_file",
                    help="File with ids from white list. If set other ids are ignored")

parser.add_argument("-m", "--max_hits_per_query", action="store", dest="max_hits_per_query", type=int, default=3,
                    help="Maximum hits per query. Default - 3")
parser.add_argument("-e", "--augustus_script_dir", action="store", dest="augustus_script_dir", default="",
                    help="Directory with augustus scripts(not set by default)")

parser.add_argument("--top_hits_priority", action="store", dest="top_hits_priority", type=int, default=20,
                    help="Priority of top exonerate hits. Default - 20")
parser.add_argument("--secondary_hits_priority", action="store", dest="secondary_hits_priority", type=int, default=10,
                    help="Priority of secondary exonerate hits. Default - 10")

parser.add_argument("--top_hits_CDS_part_cutoff", action="store", dest="top_hits_CDS_part_cutoff",
                    type=int, default=3,
                    help="CDS part cutoff for top hits, this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: 3bp ")
parser.add_argument("--secondary_hits_CDS_part_cutoff", action="store", dest="secondary_hits_CDS_part_cutoff",
                    type=int, default=9,
                    help="CDS part cutoff for secondary hits, this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: 9bp ")
parser.add_argument("--source_for_top_hits", action="store", dest="source_for_top_hits", default="EXNT",
                    help="Source for top hits. Default - 'EXNT'")
parser.add_argument("--source_for_secondary_hits", action="store", dest="source_for_secondary_hits", default="EXNS",
                    help="Source for secondary hits. Default - 'EXNS'")
parser.add_argument("-u", "--include_utr_hints", action="store_true", dest="include_utr_hints", default=False,
                    help="Include hints for UTRs. Default - False")
parser.add_argument("--max_intron_len", action="store", dest="max_intron_len", type=int,
                    help="Maximum length of intron")
parser.add_argument("--min_intron_len", action="store", dest="min_intron_len", type=int,
                    help="Minimum length of intron")

args = parser.parse_args()

top_hits_gff = "%s.target.top_hits.gff" % args.output_prefix
secondary_hits_gff = "%s.target.secondary_hits.gff" % args.output_prefix
top_hits_gff_hints = "%s.target.top_hits.hints.gff" % args.output_prefix
secondary_hits_gff_hints = "%s.target.secondary_hits.hints.gff" % args.output_prefix

Exonerate.extract_top_hits_from_target_gff(args.input, top_hits_gff, secondary_hits_gff,
                                           id_white_list_file=args.white_id_file,
                                           max_hits_per_query=args.max_hits_per_query)

AUGUSTUS.path = args.augustus_script_dir
AUGUSTUS.exonerate_to_hints(top_hits_gff, top_hits_gff_hints, priority=args.top_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.top_hits_CDS_part_cutoff, source=args.source_for_top_hits,
                            with_utrs=args.include_utr_hints)

AUGUSTUS.exonerate_to_hints(secondary_hits_gff, secondary_hits_gff_hints, priority=args.secondary_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.secondary_hits_CDS_part_cutoff, source=args.source_for_secondary_hits,
                            with_utrs=args.include_utr_hints)
