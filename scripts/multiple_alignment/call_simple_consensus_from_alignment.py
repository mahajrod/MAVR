#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.AlignInfo import SummaryInfo

from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with alignment. Default: stdin")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment file format. Default: fasta")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output")
parser.add_argument("-c", "--consensus_id", action="store", dest="consensus_id", default="consensus",
                    help="Id of consensus to use. Default: 'consensus'")
parser.add_argument("-a", "--ambiguous_symbol", action="store", dest="ambiguous_symbol", default="N",
                    help="Symbol to use for ambiguous sites. Set to X if you deal with proteins. Default: N.")
parser.add_argument("-g", "--allow_gaps", action="store_true", dest="allow_gaps", default=False,
                    help="Allow gaps in the consensus")
parser.add_argument("-e", "--cut_low_freq_regions", action="store_true", dest="cut_low_freq_regions", default=False,
                    help="Cut low frequency regions from the consensus, i.e regions where consensus ids a gap."
                         " Default: False")
parser.add_argument("-t", "--threshold", action="store", dest="threshold", default=0.5, type=float,
                    help="Threshold for the consensus. "
                         "If the most abundant nucleotide has the lower frequency, "
                         "ambiguous symbol will be used instead."
                         "Default: 0.5")

args = parser.parse_args()

alignment = AlignIO.read(FileRoutines.metaopen(args.input, "r"), format=args.format)

summary = SummaryInfo(alignment)
if args.allow_gaps or args.cut_low_freq_regions:
    consensus = summary.gap_consensus(threshold=args.threshold, ambiguous=args.ambiguous_symbol)
    if args.cut_low_freq_regions:
        consensus = consensus.replace("-", "")
else:
    consensus = summary.dumb_consensus(threshold=args.threshold, ambiguous=args.ambiguous_symbol)

SeqIO.write(SeqRecord(seq=consensus, id=args.consensus_id), args.output, format="fasta")
