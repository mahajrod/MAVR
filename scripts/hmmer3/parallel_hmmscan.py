#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.HMMER import HMMER3
from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file")
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

args = parser.parse_args()


HMMER3.threads = 1
HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output, num_of_seqs_per_scan=None, split_dir="splited_fasta",
                        splited_output_dir=args.hmmscan_output_dir, threads=args.threads,
                        combine_output_to_single_file=args.combine_output, dont_output_alignments=args.no_alignment,
                        tblout_outfile=args.tblout, domtblout_outfile=args.domtblout,
                        pfamtblout_outfile=args.pfamtblout,
                        splited_tblout_dir=args.tblout_dir, splited_domtblout_dir=args.domtblout_dir,
                        splited_pfamtblout_dir=args.pfamtblout_dir)


"""
  --acc            : prefer accessions over names in output
  --noali          : don't output alignments, so output is smaller
  --notextw        : unlimit ASCII text output line width
  --textw <n>      : set max width of ASCII text output lines  [120]  (n>=120)

Options controlling reporting thresholds:
  -E <x>     : report models <= this E-value threshold in output  [10.0]  (x>0)
  -T <x>     : report models >= this score threshold in output
  --domE <x> : report domains <= this E-value threshold in output  [10.0]  (x>0)
  --domT <x> : report domains >= this score cutoff in output

Options controlling inclusion (significance) thresholds:
  --incE <x>    : consider models <= this E-value threshold as significant
  --incT <x>    : consider models >= this score threshold as significant
  --incdomE <x> : consider domains <= this E-value threshold as significant
  --incdomT <x> : consider domains >= this score threshold as significant

Options for model-specific thresholding:
  --cut_ga : use profile's GA gathering cutoffs to set all thresholding
  --cut_nc : use profile's NC noise cutoffs to set all thresholding
  --cut_tc : use profile's TC trusted cutoffs to set all thresholding

Options controlling acceleration heuristics:
  --max    : Turn all heuristic filters off (less speed, more power)
  --F1 <x> : MSV threshold: promote hits w/ P <= F1  [0.02]
  --F2 <x> : Vit threshold: promote hits w/ P <= F2  [1e-3]
  --F3 <x> : Fwd threshold: promote hits w/ P <= F3  [1e-5]
  --nobias : turn off composition bias filter

Other expert options:
  --nonull2     : turn off biased composition score corrections
  -Z <x>        : set # of comparisons done, for E-value calculation
  --domZ <x>    : set # of significant seqs, for domain E-value calculation
  --seed <n>    : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
  --qformat <s> : assert input <seqfile> is in format <s>: no autodetection
  --daemon      : run program as a daemon
  --cpu <n>     : number of parallel CPU workers to use for multithreads
"""