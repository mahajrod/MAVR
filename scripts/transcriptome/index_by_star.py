#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Tools.Alignment import STAR
from RouToolPa.Routines.File import check_path



parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genome_dir", action="store", dest="genome_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with star index for genome")
parser.add_argument("-f", "--genome_fasta", action="store", dest="genome_fasta",
                    type=os.path.abspath,
                    help="Path to genome fasta file. If set Star will construct genome index first"
                         "in directory set by -g/--genome_dir")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use. Default - 1.")
parser.add_argument("-a", "--annotation_gtf", action="store", dest="annotation_gtf", type=os.path.abspath,
                    help="Gtf file with annotations for STAR")
parser.add_argument("-i", "--genome_size", action="store", dest="genome_size", type=int,
                    help="Genome size. Required to construct genome index")
parser.add_argument("-j", "--junction_tab_file", action="store", dest="junction_tab_file",
                    help="Junction tab file")
parser.add_argument("-r", "--star_dir", action="store", dest="star_dir", default="",
                    help="Directory with STAR binary")
parser.add_argument("-e", "--feature_from_gtf_to_use_as_exon", action="store", dest="feature_from_gtf_to_use_as_exon",
                    help="Feature from gtf to use as exon.Default: exon")
parser.add_argument("-b", "--genomeSAindexNbases", action="store", dest="genomeSAindexNbases", default=None, type=int,
                    help="Index word size. Adjust for small genomes. By default: is calculated from genome size")
parser.add_argument("-c", "--genomeChrBinNbits", action="store", dest="genomeChrBinNbits", default=None, type=int,
                    help="genomeChrBinNbits STAR option. Adjust for genomes with more than 5000 scaffolds."
                         "Use following formula min(18, log2(GenomeLength/NumberOfScaffolds))"
                         ". Default: 18")
args = parser.parse_args()

"""
Examples of usage:
skliver@dcdell:/mnt/guatemala/skliver/Pusa_sibirica/transcriptome/own/alignment/Nerpa$
~/Soft/MAVR/scripts/transcriptome/align_by_star.py \
    -d filtered/filtered/final/ \
    -g ../../genome/index_STAR \
    -t 30 \
    -o alignment/ \
    -i 2400000000 \
    -m 100000000000 \
    -r ~/Soft/STAR/bin/Linux_x86_64/

~/Soft/MAVR/scripts/transcriptome/align_by_star.py \
    -f ../../genome/nerpaGenome.1000.fa \
    -d filtered/filtered/final/ \
    -g ../../genome/index_STAR \
    -t 30 \
    -o alignment/ \
    -i 2400000000 \
    -m 100000000000 \
    -r ~/Soft/STAR/bin/Linux_x86_64/

"""

STAR.threads = args.threads
STAR.path = args.star_dir

STAR.index(args.genome_dir, args.genome_fasta, annotation_gtf=args.annotation_gtf,
           feature_from_gtf_to_use_as_exon=args.feature_from_gtf_to_use_as_exon,
           junction_tab_file=args.junction_tab_file, sjdboverhang=None,
           genomeSAindexNbases=args.genomeSAindexNbases, genomeChrBinNbits=args.genomeChrBinNbits, genome_size=args.genome_size)
