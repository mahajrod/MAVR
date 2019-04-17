#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.HMMER import HMMER3
from RouToolPa.Routines import EggNOGRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file of protein families(for example, from TreeFam)")
parser.add_argument("-s", "--input_seq", action="store", dest="input_seq",
                    help="Input file with sequences")
parser.add_argument("--no_ali", action="store_true", dest="no_alignment",
                    help="Dont save alignments to minimize output")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", default="./",
                    type=HMMER3.check_path,
                    help="Directory to write output")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-a", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for hmmer hits file. Allowed: parse, index, index_db(default)")
parser.add_argument("--hmmer_dir", action="store", dest="path", default="",
                    help="Path to directory with hmmer3.1 binaries")

args = parser.parse_args()

hits_file = "%s/%s.hits" % (args.output_dir, args.output_prefix)
top_hits_file = "%s/%s.top_hits" % (args.output_dir, args.output_prefix)
fam_file = "%s/%s.fam" % (args.output_dir, args.output_prefix)
ortholog_file = "%s/%s.orthologs" % (args.output_dir, args.output_prefix)


HMMER3.threads = 1
HMMER3.path = args.path
HMMER3.timelog = "%s/%s.timelog" % (args.output_dir, args.output_prefix)

HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output_prefix, args.output_dir,  num_of_seqs_per_scan=None,
                        threads=args.threads,
                        combine_output_to_single_file=True, dont_output_alignments=args.no_alignment,
                        biopython_165_compartibility=True
                        )

HMMER3.extract_top_hits(hits_file, args.output_prefix, parsing_mode=args.parsing_mode)
HMMER3.get_families_from_top_hits(top_hits_file, fam_file)

EggNOGRoutines.edit_profile_names_in_fam_file(fam_file, ortholog_file)
"""
for tmp_dir in args.hmmscan_output_dir, args.tblout_dir, args.domtblout_dir, args.pfamtblout_dir, "splited_fasta":
    shutil.rmtree(tmp_dir)
"""