#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Pipelines import STRPrimerPipeline

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genome_fasta", action="store", dest="genome_fasta", required=True,
                    help="Fasta file with genome sequence")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--trf_gff", action="store", dest="trf_gff",
                    help="TRF Gff file with full output including full sequence of each repeat in description field."
                         "If not set a TRF will be used to predict tandem repeats")

parser.add_argument("-k", "--directory_with_kmer_counts", action="store", dest="directory_with_kmer_counts", required=True,
                    help="Directory with files containing kmer counts")
parser.add_argument("-x", "--kmer_file_prefix", action="store", dest="kmer_file_prefix", required=True,
                    help="Prefix of files with kmer counts")
parser.add_argument("-z", "--count_kmers", action="store_true", dest="count_kmers",
                    help="Count kmers by Glistmaker for primer 3 using path specified "
                         "by -k/--directory_with_kmer_counts and -x/--kmer_file_prefix to store the results")

parser.add_argument("-n", "--min_period", action="store", dest="min_period", type=int, default=3,
                    help="Minimum period of repeats to extract. Default: 3")
parser.add_argument("-m", "--max_period", action="store", dest="max_period", type=int, default=5,
                    help="Maximum period of repeats to extract. Default: 5")
parser.add_argument("-b", "--min_copy_number", action="store", dest="min_copy_number", type=float, default=20.0,
                    help="Minimum number of copies(float) to extract. Default: 20.0")
parser.add_argument("-c", "--min_perfect_copy_number", action="store", dest="min_perfect_copy_number", type=int,
                    help="Minimum number of perfect copies to extract. Default: not set")
parser.add_argument("-j", "--require_tandem_perfect_copies", action="store_true", default=False,
                    dest="require_tandem_perfect_copies",
                    help="Require perfect copies to be tandem inside STR")

parser.add_argument("-a", "--max_copy_number", action="store", dest="max_copy_number", type=float,
                    help="Maximum number of copies(float) to extract. Default: not set")
parser.add_argument("-p", "--pattern", action="store", dest="pattern",
                    help="Extract patterns only with this pattern")
parser.add_argument("-d", "--min_percentage_of_matches", action="store", dest="min_percentage_of_matches", type=int,
                    help="Minimum percentage of matches(int). Default: not set")
parser.add_argument("-e", "--max_percentage_of_indels", action="store", dest="max_percentage_of_indels", type=int,
                    help="Maximum percentage of indels(int). Default: not set")

parser.add_argument("-l", "--left_flank_len", action="store", dest="left_flank_len", type=int, default=200,
                    help="Length of left flank. Default: 200")
parser.add_argument("-r", "--right_right_len", action="store", dest="right_flank_len", type=int, default=200,
                    help="Length of right flank. Default: 200")

parser.add_argument("--min_gap_len", action="store", dest="min_gap_len", default=5, type=int,
                    help="Minimum length of polyN to be treated as gap. Default: 5")

parser.add_argument("--primer3_dir", action="store", dest="primer3_dir", default="",
                    help="Directory with primer3_core binary")
parser.add_argument("--glistmaker_dir", action="store", dest="glistmaker_dir", default="",
                    help="Directory with Glistmaker binary")
parser.add_argument("--trf_dir", action="store", dest="trf_dir", default="",
                    help="Directory with TRF binary")

parser.add_argument("--primer3_thermo_config_dir", action="store", dest="primer3_thermo_config_dir",
                    help="Directory with primer3 config for thermodynamic approach")

parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default: 1")
"""
parser.add_argument("-s", "--split_output_by_monomer_len", action="store_true",
                    dest="split_output_by_monomer_len", default=False,
                    help="Split output by STR monomer length")
"""
args = parser.parse_args()

"""
#usage example(with only genome sequence as input, so trf and glistmaker will be started):
time ~/Soft/MAVR/scripts/primers/generate_str_primers.py -o assembly.hybrid.all.trf \
                                                          -g ../../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta \
                                                          -k ./kmers/ \
                                                          -x mustela_nigripes_hybrid \
                                                          --primer3_dir ~/Soft/primer3/src/ \
                                                          --primer3_thermo_config_dir ~/Soft/primer3/src/primer3_config/ \
                                                          -t 30 \
                                                          --glistmaker_dir ~/Soft/GenomeTester4/bin/ \
                                                          --trf_dir ~/Soft/TRF/ \
                                                          -c 20 -b 20 -j -n 3 -m 5 -z \

#usage with trf gff as input with precounted kmer frequences:
time ~/Soft/MAVR/scripts/primers/generate_str_primers.py -o assembly.hybrid.all.trf \
                                                         -f ../assembly.hybrid.all.trf.with_rep_seqs.gff \
                                                         -g ../../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta \
                                                         -k ../../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ \
                                                         -x mustela_nigripes_hybrid \
                                                         --primer3_dir ~/Soft/primer3/src/ \
                                                         --primer3_thermo_config_dir ~/Soft/primer3/src/primer3_config/ \
                                                         -t 30 \
                                                         --glistmaker_dir ~/Soft/GenomeTester4/bin/ \
                                                         --trf_dir ~/Soft/TRF/ \
                                                         -c 20 -b 20 -j -n 3 -m 5

real	1m47.007s
user	1m26.900s
sys	0m25.032s

time ~/Soft/MAVR/scripts/primers/generate_str_primers.py -o assembly.hybrid.all.trf \
                                                         -f ../assembly.hybrid.all.trf.with_rep_seqs.gff \
                                                         -g ../../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta \
                                                         -k ../../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ \
                                                         -x mustela_nigripes_hybrid \
                                                         --primer3_dir ~/Soft/primer3/src/ \
                                                         --primer3_thermo_config_dir ~/Soft/primer3/src/primer3_config/ \
                                                         -t 30 \
                                                         --glistmaker_dir ~/Soft/GenomeTester4/bin/ \
                                                         --trf_dir ~/Soft/TRF/ \
                                                         -c 12 -b 12 -j -n 3 -m 6

"""

STRPrimerPipeline.threads = args.threads

STRPrimerPipeline.trf_dir = args.trf_dir
STRPrimerPipeline.primer3_dir = args.primer3_dir
STRPrimerPipeline.glistmaker_dir = args.glistmaker_dir

STRPrimerPipeline.primer3_thermo_config_dir = args.primer3_thermo_config_dir

STRPrimerPipeline.primer_prediction_pipeline(args.genome_fasta, args.output_prefix, trf_gff=args.trf_gff,
                                             min_str_period=args.min_period, max_str_period=args.max_period,
                                             min_copy_number=args.min_copy_number,
                                             max_copy_number=args.max_copy_number, pattern=None,
                                             min_perfect_copy_number=args.min_perfect_copy_number,
                                             require_tandem_perfect_copies=args.require_tandem_perfect_copies,
                                             left_flank_len=args.left_flank_len, right_flank_len=args.right_flank_len,
                                             core_seq_coords_entry="core_seq_coords",
                                             id_description_entry="ID",
                                             kmer_dir=args.directory_with_kmer_counts,
                                             kmer_file_prefix=args.kmer_file_prefix,
                                             count_kmers=args.count_kmers,
                                             min_percentage_of_matches=args.min_percentage_of_matches,
                                             max_percentage_of_indels=args.max_percentage_of_indels,
                                             optimal_primer_len=None,
                                             min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                                             softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                                             optimal_melting_temperature=None, min_melting_temperature=None,
                                             max_melting_temperature=None, black_list_of_seqs_fasta=None,
                                             trf_matching_weight=2, trf_mismatching_penalty=7,
                                             trf_indel_penalty=7, trf_matching_probability=80, trf_indel_probability=10,
                                             trf_min_score=50, trf_max_period_size=500, threads=None,
                                             min_gap_len=args.min_gap_len)

