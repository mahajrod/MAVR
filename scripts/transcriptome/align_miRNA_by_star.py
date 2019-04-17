#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Tools.Alignment import STAR
from RouToolPa.Tools.Samtools import SamtoolsV1
from Pipelines import Pipeline
from RouToolPa.Routines import FileRoutines
from RouToolPa.Routines.File import check_path



parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-g", "--genome_dir", action="store", dest="genome_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with star index for genome")
parser.add_argument("-f", "--genome_fasta", action="store", dest="genome_fasta",
                    type=os.path.abspath,
                    help="Path to genome fasta file. If set Star will construct genome index first"
                         "in directory set by -g/--genome_dir")
parser.add_argument("-s", "--samples", action="store", dest="samples", type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use in Trimmomatic. Default - 1.")
parser.add_argument("-a", "--annotation_gtf", action="store", dest="annotation_gtf", type=os.path.abspath,
                    help="Gtf file with annotations for STAR")
parser.add_argument("-i", "--genome_size", action="store", dest="genome_size", type=int,
                    help="Genome size. Required for constructing genome index")

parser.add_argument("-r", "--star_dir", action="store", dest="star_dir", default="",
                    help="Directory with STAR binary")

parser.add_argument("-m", "--max_memory_for_bam_sorting", action="store", type=int,
                    dest="max_memory_for_bam_sorting", default=8000000000,
                    help="Max memory for bam sorting. Default: 8 000 000 000")
parser.add_argument("-e", "--enable_soft_clipping", action="store_false", default=True,
                    dest="enable_soft_clipping",
                    help="Enable soft clipping. Default: False")
parser.add_argument("-x", "--max_number_of_alignments_per_read", action="store", type=int,
                    dest="max_number_of_alignments_per_read", default=10,
                    help="Maximum number of alignments per read. Default: 10")
parser.add_argument("-y", "--max_number_of_mismatches", action="store", type=int,
                    dest="max_number_of_mismatches", default=1,
                    help="Maximum number of mismatches. Ignored if set -z/--max_relative number_of_mismatches. "
                         "Default: 1")
parser.add_argument("-z", "--max_relative_number_of_mismatches", action="store", type=float,
                    dest="max_relative_number_of_mismatches",
                    help="Maximum number of mismatches relative to read length. "
                         "If set overrides -y/--max_number_of_mismatches option. Default: Not set")
max_relative_number_of_mismatches=None

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

if args.genome_fasta:
    STAR.index(args.genome_dir, args.genome_fasta, annotation_gtf=args.annotation_gtf,
               junction_tab_file=args.junction_tab_file, sjdboverhang=None,
               genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=args.genome_size)

sample_list = args.samples if args.samples else Pipeline.get_sample_list(args.samples_dir)

FileRoutines.safe_mkdir(args.output_dir)

for sample in sample_list:
    print ("Handling %s" % sample)
    sample_dir = "%s/%s/" % (args.samples_dir, sample)
    alignment_sample_dir = "%s/%s/" % (args.output_dir, sample)
    FileRoutines.safe_mkdir(alignment_sample_dir)
    filetypes, forward_files, reverse_files, se_files = FileRoutines.make_lists_forward_and_reverse_files(sample_dir)

    print("\tAligning reads...")

    STAR.align_miRNA(args.genome_dir, se_files,
                     output_dir=alignment_sample_dir,
                     annotation_gtf=args.annotation_gtf if not args.genome_fasta else None,
                     max_memory_for_bam_sorting=args.max_memory_for_bam_sorting,
                     max_alignments_per_read=args.max_number_of_alignments_per_read,
                     no_soft_clip=args.enable_soft_clipping,
                     max_number_of_mismatches=args.max_number_of_mismatches,
                     max_relative_number_of_mismatches=args.max_relative_number_of_mismatches)

    print("\tIndexing bam file...")
    resulting_bam_file = "%s/Aligned.sortedByCoord.out.bam" % alignment_sample_dir
    SamtoolsV1.index(resulting_bam_file)