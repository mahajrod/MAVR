#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.HMMER import HMMER3
from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file of protein families(for example, from TreeFam)")
parser.add_argument("-s", "--input_seq", action="store", dest="input_seq",
                    help="Input file with sequences")
parser.add_argument("--no_ali", action="store_true", dest="no_alignment",
                    help="Dont save alignments to minimize output")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-d", "--hmmscan_output_dir", action="store", dest="hmmscan_output_dir",
                    default="hmmscan_output_dir/", type=check_path,
                    help="Directory to write intermediate(splited) output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("--tblout_dir", action="store", dest="tblout_dir",
                    default="tblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-sequence hits")
parser.add_argument("--domtblout_dir", action="store", dest="domtblout_dir",
                    default="domtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-domain hits")
parser.add_argument("--pfamtblout_dir", action="store", dest="pfamtblout_dir",
                    default="pfamtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) table of hits and domains to file, in Pfam format ")
parser.add_argument("--hmmer_dir", action="store", dest="path", default="",
                    help="Path to directory with hmmer3.1 binaries")

args = parser.parse_args()

output_dir = "./"
hits_file = "%s.hits" % args.output_prefix
top_hits_file = "%s.top_hits" % args.output_prefix
not_significant_ids_file = "%s.not_significant" % args.output_prefix
not_found_ids_file = "%s.not_found" % args.output_prefix
fam_file = "%s.fam" % args.output_prefix

HMMER3.threads = 1
HMMER3.path = args.path
HMMER3.timelog = "%s.timelog" % args.output_prefix

HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output_prefix, output_dir, num_of_seqs_per_scan=None, split_dir="splited_fasta",
                        splited_output_dir=args.hmmscan_output_dir, threads=args.threads,
                        combine_output_to_single_file=True, dont_output_alignments=args.no_alignment,
                        splited_tblout_dir=args.tblout_dir, splited_domtblout_dir=args.domtblout_dir,
                        splited_pfamtblout_dir=args.pfamtblout_dir,
                        biopython_165_compartibility=True
                        )

HMMER3.extract_top_hits(hits_file, top_hits_file, not_significant_ids_file=not_significant_ids_file,
                        not_found_ids_file=not_found_ids_file)
HMMER3.get_families_from_top_hits(top_hits_file, fam_file)
