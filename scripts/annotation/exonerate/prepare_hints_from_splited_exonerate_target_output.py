#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_prefix", action="store", dest="input_prefix", required=True,
                    help="Prefix of files with splited exonerate target output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-d", "--white_id_file", action="store", dest="white_id_file",
                    help="File with ids from white list. If set other ids are ignored")
"""
parser.add_argument("-m", "--max_hits_per_query", action="store", dest="max_hits_per_query", type=int, default=3,
                    help="Maximum hits per query. Default - 3")
"""

parser.add_argument("-e", "--augustus_script_dir", action="store", dest="augustus_script_dir", default="",
                    help="Directory with augustus scripts(not set by default)")

parser.add_argument("--other_top_hits_priority", action="store", dest="other_top_hits_priority", type=int, default=20,
                    help="Priority of top(not full) exonerate hits. Default - 20")
parser.add_argument("--other_secondary_hits_priority", action="store", dest="other_secondary_hits_priority", type=int, default=10,
                    help="Priority of secondary(not full) exonerate hits. Default - 10")
parser.add_argument("--full_top_hits_priority", action="store", dest="full_top_hits_priority", type=int, default=100,
                    help="Priority of top(full) exonerate hits. Default - 100")
parser.add_argument("--full_secondary_hits_priority", action="store", dest="full_secondary_hits_priority", type=int, default=10,
                    help="Priority of secondary(full) exonerate hits. Default - 10")

parser.add_argument("--other_top_hits_CDS_part_cutoff", action="store", dest="other_top_hits_CDS_part_cutoff",
                    type=int, default=3,
                    help="CDS part cutoff for top hits(not full), this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: 3bp ")
parser.add_argument("--other_secondary_hits_CDS_part_cutoff", action="store", dest="other_secondary_hits_CDS_part_cutoff",
                    type=int, default=9,
                    help="CDS part cutoff for secondary(not full) hits, this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: 9bp ")
parser.add_argument("--full_top_hits_CDS_part_cutoff", action="store", dest="full_top_hits_CDS_part_cutoff",
                    type=int, default=0,
                    help="CDS part cutoff for top hits(full), this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: no cutoff ")
parser.add_argument("--full_secondary_hits_CDS_part_cutoff", action="store", dest="full_secondary_hits_CDS_part_cutoff",
                    type=int, default=3,
                    help="CDS part cutoff for secondary(full) hits, this many bp are cut off of each CDSpart hint w.r.t. "
                         "the exonerate cds. Default: 3bp ")

parser.add_argument("--source_for_other_top_hits", action="store", dest="source_for_other_top_hits", default="EXNT",
                    help="Source for top(not full) hits. Default - 'EXNT'")
parser.add_argument("--source_for_other_secondary_hits", action="store", dest="source_for_other_secondary_hits", default="EXNS",
                    help="Source for secondary(not full) hits. Default - 'EXNS'")
parser.add_argument("--source_for_full_top_hits", action="store", dest="source_for_full_top_hits", default="EXNFULLT",
                    help="Source for top(full) hits. Default - 'EXNT'")
parser.add_argument("--source_for_full_secondary_hits", action="store", dest="source_for_full_secondary_hits", default="EXNS",
                    help="Source for secondary(full) hits. Default - 'EXNS'")

parser.add_argument("-u", "--include_utr_hints", action="store_true", dest="include_utr_hints", default=False,
                    help="Include hints for UTRs. Default - False")
parser.add_argument("--max_intron_len", action="store", dest="max_intron_len", type=int,
                    help="Maximum length of intron")
parser.add_argument("--min_intron_len", action="store", dest="min_intron_len", type=int,
                    help="Minimum length of intron")

args = parser.parse_args()

"""
Example:
~/Soft/MAVR/scripts/annotation/prepare_hints_from_splited_exonerate_target_output.py -i mustela_putoris_furo \
                                                                                     -e ~/Soft/augustus-3.2.1/scripts/ \
                                                                                     -o mustela_putoris_furo
"""

input_precise_top_gff = "%s.precise_top.target.gff" % args.input_prefix
input_precise_secondary_gff = "%s.precise_secondary.target.gff" % args.input_prefix

input_other_top_gff = "%s.other_top.target.gff" % args.input_prefix
input_other_secondary_gff = "%s.other_secondary.target.gff" % args.input_prefix

precise_top_hints = "%s.precise_top.target.hints.gff" % args.output_prefix
precise_secondary_hints = "%s.precise_secondary.target.hints.gff" % args.output_prefix

other_top_hints = "%s.other_top.target.hints.gff" % args.output_prefix
other_secondary_hints = "%s.other_secondary.target.hints.gff" % args.output_prefix


AUGUSTUS.path = args.augustus_script_dir

AUGUSTUS.exonerate_to_hints(input_precise_top_gff, precise_top_hints, priority=args.full_top_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.full_top_hits_CDS_part_cutoff, source=args.source_for_full_top_hits,
                            with_utrs=args.include_utr_hints)

AUGUSTUS.exonerate_to_hints(input_precise_secondary_gff, precise_secondary_hints,
                            priority=args.full_secondary_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.full_secondary_hits_CDS_part_cutoff,
                            source=args.source_for_full_secondary_hits,
                            with_utrs=args.include_utr_hints)

AUGUSTUS.exonerate_to_hints(input_other_top_gff, other_top_hints, priority=args.other_top_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.other_top_hits_CDS_part_cutoff, source=args.source_for_other_top_hits,
                            with_utrs=args.include_utr_hints)

AUGUSTUS.exonerate_to_hints(input_other_secondary_gff, other_secondary_hints,
                            priority=args.other_secondary_hits_priority,
                            min_intron_len=args.min_intron_len, max_intron_len=args.max_intron_len,
                            CDS_part_cutoff=args.other_secondary_hits_CDS_part_cutoff,
                            source=args.source_for_other_secondary_hits,
                            with_utrs=args.include_utr_hints)
