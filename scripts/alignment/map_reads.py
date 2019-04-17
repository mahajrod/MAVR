#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Samtools import SamtoolsV1
from RouToolPa.Tools.Picard import AddOrReplaceReadGroups, MarkDuplicates


def make_list_from_comma_sep_string(s):
    return s.split(",")

parser = argparse.ArgumentParser()

parser.add_argument("-p", "--prefix", action="store", dest="prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-n", "--sample_name", action="store", dest="sample_name", required=True,
                    help="Sample name")
parser.add_argument("-r", "--forward_reads", action="store", dest="forward_reads",
                    type=make_list_from_comma_sep_string,
                    help="Comma-separated list of files with forward reads")
parser.add_argument("-l", "--reverse_reads", action="store", dest="reverse_reads",
                    type=make_list_from_comma_sep_string,
                    help="Comma-separated list of files with reverse reads")
parser.add_argument("-u", "--unpaired_reads", action="store", dest="unpaired_reads",
                    type=make_list_from_comma_sep_string,
                    help="Comma-separated list of files with unpaired reads")

parser.add_argument("-i", "--aligner_index", action="store", dest="index",
                    help="Aligner-specific index")

parser.add_argument("-a", "--aligner", action="store", dest="aligner", default="bwa",
                    help="Aligner to use. Possible aligners: bwa(default), bowtie2")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")

parser.add_argument("-q", "--quality", action="store", dest="quality", default="phred33",
                    help="Quality type. Possible variants phred33, phred64. Default: phred33")
parser.add_argument("-j", "--alignment_format", action="store", dest="alignment_format", default="bam",
                    help="Format of output alignments. Allowed: bam(default), sam, cram")


parser.add_argument("-g", "--add_read_groups_by_picard", action="store_true", dest="add_read_groups_by_picard", default=False,
                    help="Add read groups to final bam using PICARD. Use this option if aligner don't support adding readgroups itself")
parser.add_argument("-d", "--picard_dir", action="store", dest="picard_dir", default="",
                    help="Path to Picard directory. Required to add read groups")
parser.add_argument("-e", "--retain_intermediate_files", action="store_true", dest="retain_temp", default=False,
                    help="Retain intermediate files")
"""
parser.add_argument("-z", "--calculate_median_coverage", action="store_true", dest="calculate_median_coverage",
                    default=False,
                    help="Calculate median coverage")
parser.add_argument("-x", "--calculate_mean_coverage", action="store_true", dest="calculate_mean_coverage",
                    default=False,
                    help="Calculate mean coverage")
parser.add_argument("-f", "--flanks_size", action="store", dest="flanks_size",
                    default=0, type=int,
                    help="Size of flanks to remove when calculating mean/median coverage")
"""



parser.add_argument("-b", "--bed_with_regions", action="store", dest="bed",
                    help="Bed file with regions to output.")

parser.add_argument("-c", "--black_flag_value", action="store", dest="black_flag_value", type=int,
                    help="Black bam flag value. By default unaligned, supplementary alignments and "
                         "nonprimary alignments will be removed")

parser.add_argument("-e", "--white_flag_value", action="store", dest="white_flag_value", type=int,
                    help="White flag value")
parser.add_argument("-m", "--max_memory", action="store", dest="max_memory", default="30g",
                    help="Maximum memory to use. Default: 30g")

args = parser.parse_args()

black_flag_value = args.black_flag_value if args.black_flag_value else \
    SamtoolsV1.bam_flags["read_unmapped"] + SamtoolsV1.bam_flags["supplementary_alignment"] \
    + SamtoolsV1.bam_flags["not_primary_alignment"]

sorted_alignment = "%s.%s" % (args.prefix, args.alignment_format)
sorted_alignment_picard_groups = None

final_alignment = "%s.mkdup.%s" % (args.prefix, args.alignment_format)
duplicates_stat_file = "%s.duplicates.stat" % args.prefix
coverage_file = "%s.coverage.bed" % args.prefix

if args.aligner == "bowtie2":
    aligner = Bowtie2
elif args.aligner == "bwa":
    aligner = BWA
    
else:
    raise ValueError("")

aligner.threads = args.threads
MarkDuplicates.jar_path = args.picard_dir
MarkDuplicates.max_memory = args.max_memory
AddOrReplaceReadGroups.jar_path = args.picard_dir

aligner.align(args.index, forward_reads_list=args.forward_reads, reverse_reads_list=args.reverse_reads,
              unpaired_reads_list=args.unpaired_reads, quality_score=args.quality, output_prefix=args.prefix,
              output_format=args.alignment_format,
              read_group_name=args.sample_name,
              PU="x",
              SM=args.sample_name,
              platform="Illumina",
              LB="x",
              sort_by_coordinate=True,
              sort_by_name=False,
              max_per_sorting_thread_memory="10G")


if args.add_read_groups_by_picard:
    sorted_alignment_picard_groups = "%s.picard_groups.%s" % (args.prefix, args.alignment_format)
    AddOrReplaceReadGroups.add_read_groups(sorted_alignment, sorted_alignment_picard_groups,
                                           RGID=args.prefix, RGLB=args.prefix, RGPL=args.prefix,
                                           RGSM=args.prefix, RGPU=args.prefix)

if args.alignment_format == "bam":
    SamtoolsV1.index(sorted_alignment_picard_groups if sorted_alignment_picard_groups else sorted_alignment)

MarkDuplicates.run(sorted_alignment_picard_groups if sorted_alignment_picard_groups else sorted_alignment,
                   final_alignment,
                   duplicates_stat_file)

if args.alignment_format == "bam":
    SamtoolsV1.index(final_alignment)

"""
GenomeCov.get_coverage(final_alignment, genome_bed, coverage_file)
if not args.retain_temp:
    os.remove(sorted_alignment)
    if args.add_read_groups_by_picard:
        os.remove(sorted_alignment_picard_groups)

if args.calculate_median_coverage or args.calculate_mean_coverage:
    coverage_dict = SynDict()
    coverage_dict.read(coverage_file, header=False, separator="\t", allow_repeats_of_key=True,
                       values_separator=",", key_index=0, value_index=2, expression=int)
    if args.calculate_median_coverage:
        with open("%s_median_coverage.tab" % args.prefix, "w") as out_fd:
            for region in coverage_dict:
                mediana = median(array(coverage_dict[region] if args.flanks_size == 0
                                       else coverage_dict[region][args.flanks_size:-args.flanks_size]))
                out_fd.write("%s\t%f\n" % (region, mediana))
    if args.calculate_mean_coverage:
        with open("%s_mean_coverage.tab" % args.prefix, "w") as out_fd:
            for region in coverage_dict:
                meana = mean(array(coverage_dict[region] if args.flanks_size == 0
                                   else coverage_dict[region][args.flanks_size:-args.flanks_size]))
                out_fd.write("%s\t%f\n" % (region, meana))
"""