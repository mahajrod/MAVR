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
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file")
parser.add_argument("-c", "--combine_output", action="store_true", dest="combine_output",
                    help="Combine output files to single")
parser.add_argument("--no_ali", action="store_true", dest="no_alignment",
                    help="Dont save alignments to minimize output")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-d", "--hmmscan_output_dir", action="store", dest="hmmscan_output_dir",
                    default="hmmscan_output_dir/", type=check_path,
                    help="Directory to write intermediate(splited) output")
parser.add_argument("-p", "--top_hits_file", action="store", dest="top_hits_file",
                    default="hits.top_hits",
                    help="File to write top hits")
parser.add_argument("-n", "--not_found_ids_file", action="store", dest="not_found_file",
                    default="not_found.ids",
                    help="Ids of not found proteins")
parser.add_argument("-a", "--not_significant_ids_file", action="store", dest="not_significant_file",
                    default="not_significant.ids",
                    help="Ids of proteins with not_significant hits")
parser.add_argument("-f", "--fam_file", action="store", dest="fam_file",
                    default="proteins.fam",
                    help="File to write families of proteins")

parser.add_argument("--tblout_dir", action="store", dest="tblout_dir",
                    default="tblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-sequence hits")
parser.add_argument("--domtblout_dir", action="store", dest="domtblout_dir",
                    default="domtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-domain hits")
parser.add_argument("--pfamtblout_dir", action="store", dest="pfamtblout_dir",
                    default="pfamtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) table of hits and domains to file, in Pfam format ")
parser.add_argument("--tblout", action="store", dest="tblout",
                    help="File to save parseable table of per-sequence hits")
parser.add_argument("--domtblout", action="store", dest="domtblout",
                    help="File to save parseable table of per-domain hits")
parser.add_argument("--pfamtblout", action="store", dest="pfamtblout",
                    help="File to save table of hits and domains to file, in Pfam format ")
parser.add_argument("--hmmer_dir", action="store", dest="path", default="",
                    help="Path to directory with hmmer3.1 binaries")

args = parser.parse_args()


HMMER3.threads = 1
HMMER3.path = args.path

""""
HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output, num_of_seqs_per_scan=None, split_dir="splited_fasta",
                        splited_output_dir=args.hmmscan_output_dir, threads=args.threads,
                        combine_output_to_single_file=args.combine_output, dont_output_alignments=args.no_alignment,
                        tblout_outfile=args.tblout, domtblout_outfile=args.domtblout,
                        pfamtblout_outfile=args.pfamtblout,
                        splited_tblout_dir=args.tblout_dir, splited_domtblout_dir=args.domtblout_dir,
                        splited_pfamtblout_dir=args.pfamtblout_dir
                        )
"""
HMMER3.extract_top_hits(args.output, args.top_hits_file, not_significant_ids_file=args.not_significant_file,
                        not_found_ids_file=args.not_found_file)
HMMER3.get_families_from_top_hits(args.top_hits_file, args.fam_file)
