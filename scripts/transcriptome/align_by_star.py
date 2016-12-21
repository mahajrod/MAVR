#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Tools.Alignment import STAR
from Pipelines import Pipeline
from Routines import FileRoutines
from Routines.File import check_path

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
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use in Trimmomatic. Default - 1.")

parser.add_argument("-a", "--annotation_gtf", action="store", dest="annotation_gtf", type=os.path.abspath,
                    help="Gtf file with annotations for STAR")

parser.add_argument("-i", "--genome_size", action="store", dest="genome_size", type=int,
                    help="Genome size. Required for constructing genome index")
parser.add_argument("-j", "--junction_tab_file", action="store", dest="junction_tab_file",
                    help="Junction tab file")
parser.add_argument("-r", "--star_dir", action="store", dest="star_dir", default="",
                    help="Directory with STAR binary")
parser.add_argument("-u", "--include_unmapped_reads", action="store_true",
                    dest="include_unmapped_reads", default=False,
                    help="Include unmapped reads in Bam file")
parser.add_argument("-m", "--max_memory_for_bam_sorting", action="store", type=int,
                    dest="max_memory_for_bam_sorting", default=8000000000,
                    help="Max memory for bam sorting")
parser.add_argument("-x", "--max_intron_length", action="store", dest="max_intron_length", type=int,
                    help="Maximum intron length. Default: not set")

args = parser.parse_args()

STAR.threads = args.threads
STAR.path = args.star_dir

if args.genome_fasta:
    STAR.index(args.genome_dir, args.genome_fasta, annotation_gtf=args.annotation_gtf,
               junction_tab_file=args.junction_tab_file, sjdboverhang=None,
               genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=args.genome_size)

sample_list = args.samples if args.samples else Pipeline.get_sample_list(args.samples_directory)


for sample in sample_list:
    print ("Handling %s" % sample)
    sample_dir = "%s/%s/" % (args.samples_directory, sample)
    alignment_sample_dir = "%s/%s/" % (args.output_dir, sample)
    filetypes, forward_files, reverse_files = FileRoutines.make_lists_forward_and_reverse_files(sample_dir)

    print "\tAligning reads..."

    STAR.align(args.genome_dir, forward_files, reverse_read_list=reverse_files, annotation_gtf=annotation_gtf,
               feature_from_gtf_to_use_as_exon=feature_from_gtf_to_use_as_exon,
               exon_tag_to_use_as_transcript_id=exon_tag_to_use_as_transcript_id,
               exon_tag_to_use_as_gene_id=exon_tag_to_use_as_gene_id,
               length_of_sequences_flanking_junction=length_of_sequences_flanking_junction,
               junction_tab_file_list=junction_tab_file_list,
               three_prime_trim=three_prime_trim, five_prime_trim=five_prime_trim,
               adapter_seq_for_three_prime_clip=adapter_seq_for_three_prime_clip,
               max_mismatch_percent_for_adapter_trimming=max_mismatch_percent_for_adapter_trimming,
               three_prime_trim_after_adapter_clip=three_prime_trim_after_adapter_clip,
               output_type=output_type, sort_bam=sort_bam,
               max_memory_for_bam_sorting=max_memory_for_bam_sorting,
               include_unmapped_reads_in_bam=args.include_unmapped_reads,
               output_unmapped_reads=output_unmapped_reads, output_dir=alignment_sample_dir,
               two_pass_mode=two_pass_mode, max_intron_length=max_intron_length)