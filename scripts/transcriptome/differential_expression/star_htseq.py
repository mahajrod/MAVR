#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Pipelines import DiffExpressionPipeline
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
parser.add_argument("-e", "--gff_for_htseq", action="store", required=True,
                    dest="gff_for_htseq", help="Gff file with annotations for HTseq")
parser.add_argument("-b", "--count_table_file_prefix", action="store", dest="count_table_file_prefix",
                    type=os.path.abspath,
                    help="Prefix of files to write resulting count table")
parser.add_argument("-n", "--stranded_rnaseq", action="store", default="yes",
                    dest="stranded_rnaseq", help="Type of RNAseq data. Allowed: 'yes' - stranded"
                                                 "'no' - unstranded, 'reverse' - stranded but with "
                                                 "reversed orientation of reads in pair. Default - yes.")
parser.add_argument("-l", "--min_alignment_quality", action="store", default=10, type=int,
                    dest="min_alignment_quality",
                    help="Minimum quality of read alignment to be considered. Default - 10")

parser.add_argument("--feature_type_for_htseq", action="store", default="exon",
                    dest="feature_type_for_htseq",
                    help="Feature type from annotation gff to be used for counting reads. Default - 'exon'")
parser.add_argument("--feature_id_attribute_for_htseq", action="store", default="gene_id",
                    dest="feature_id_attribute_for_htseq",
                    help="Feature id attribute from annotation gff to be considered for counting reads. "
                         "Default - 'gene_id'")
parser.add_argument("--htseq_mode", action="store", default="union",
                    dest="htseq_mode",
                    help="HTSeq mode for counting reads. Default - 'union'")
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
                    help="Gtf file with annotations")
parser.add_argument("-i", "--genome_size", action="store", dest="genome_size", type=int,
                    help="Genome size. Required for constructing genome index")
parser.add_argument("-j", "--junction_tab_file", action="store", dest="junction_tab_file",
                    help="Junction tab file")
parser.add_argument("-r", "--star_dir", action="store", dest="star_dir", default="",
                    help="Directory with STAR binary")
parser.add_argument("-u", "--include_unmapped_reads", action="store_true",
                    dest="include_unmapped_reads", default=False,
                    help="Include unmapped reads in Bam file")
parser.add_argument("-m", "--max_memory_per_thread_for_bam_sorting", action="store",
                    dest="max_memory_per_thread_for_bam_sorting", default="4G",
                    help="Max memory per thread for bam sorting. Default: 4G")
parser.add_argument("-x", "--max_intron_length", action="store", dest="max_intron_length", type=int,
                    help="Maximum intron length. Default: not set")

""""

parser.add_argument("-k", "--adapter_kmers", action="store", dest="adapter_kmers", type=os.path.abspath,
                    required=True,
                    help="File with adapter k-mers for Coockiecutter")

parser.add_argument("-p", "--pe_score", action="store", dest="pe_score", type=int, default=30,
                    help="PE reads adapter score. Works only if -a/--adapters option is set. Default - 30.")
parser.add_argument("-e", "--se_score", action="store", dest="se_score", type=int, default=10,
                    help="SE reads adapter score. Works only if -a/--adapters option is set. Default - 10.")
parser.add_argument("-n", "--min_adapter_len", action="store", dest="min_adapter_len", type=int, default=1,
                    help="Minimum length of adapter fragment. Works only if -a/--adapters option is set. Default - 1.")

parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=20,
                    type=int,
                    help="Quality threshold for sliding window or whole read."
                         "Depends on -q/--average_quality_threshold option.Default - 15.")
parser.add_argument("-b", "--base_quality", action="store", dest="base_quality", default="phred33",
                    help="Type of base quality. Possible variants: phred33, phred64. Default - phred33 ")

parser.add_argument("-l", "--min_length", action="store", dest="min_len", type=int, default=50,
                    help="Minimum length of read to retain. Default - 50")


parser.add_argument("-c", "--coockiecutter_dir", action="store", dest="coockiecutter_dir", default="",
                    help="Path to Coockiecutter directory")

"""
args = parser.parse_args()

"""
EXAMPLE
cd ~/workdir/yeast/nizhnikov/good_run
~/soft/MAVR/scripts/transcriptome/differential_expression/star_htseq.py -d fastq/filtered/final/ -o ./ \
                                                                        -g ~/data/genomes/saccharomyces_cerevisiae/S288C_R64/fasta/STAR_index \
                                                                        -f ~/data/genomes/saccharomyces_cerevisiae/S288C_R64/fasta/S288C_reference_sequence_R64-1-1_20110203_modified.fasta \
                                                                        -t 10 -i 12000000 \
                                                                        -j ~/data/genomes/saccharomyces_cerevisiae/S288C_R64/gff/SJ.intron.tab -r ~/soft/STAR/bin/Linux_x86_64// \
                                                                        -m 20000000000 \
                                                                        -x 3000 \
                                                                        -n no \
                                                                        -e ~/data/genomes/saccharomyces_cerevisiae/S288C_R64/gff/saccharomyces_cerevisiae_R64-1-1_20110208_edited_mt_no_fasta.gff \
                                                                        --feature_type_for_htseq CDS \
                                                                        --feature_id_attribute_for_htseq Parent \
                                                                        -b all_samples_read_count.table
"""
DiffExpressionPipeline.star_and_htseq(args.genome_dir, args.samples_dir, args.output_dir, args.gff_for_htseq,
                                      args.count_table_file_prefix, genome_fasta=args.genome_fasta,
                                      genome_size=args.genome_size, samples_to_handle=args.samples,
                                      annotation_gtf=args.annotation_gtf,
                                      feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                                      exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                                      junction_tab_file_list=args.junction_tab_file,
                                      three_prime_trim=None, five_prime_trim=None,
                                      adapter_seq_for_three_prime_clip=None,
                                      max_mismatch_percent_for_adapter_trimming=None,
                                      three_prime_trim_after_adapter_clip=None,
                                      output_type="BAM", sort_bam=True,
                                      max_memory_per_thread_for_bam_sorting=args.max_memory_per_thread_for_bam_sorting,
                                      include_unmapped_reads_in_bam=args.include_unmapped_reads,
                                      output_unmapped_reads=True,  two_pass_mode=False, star_dir=args.star_dir,
                                      threads=args.threads, max_intron_length=args.max_intron_length,
                                      stranded_rnaseq=args.stranded_rnaseq,
                                      min_alignment_quality=args.min_alignment_quality,
                                      feature_type_for_htseq=args.feature_type_for_htseq,
                                      feature_id_attribute_for_htseq=args.feature_id_attribute_for_htseq,
                                      htseq_mode=args.htseq_mode)
