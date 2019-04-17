#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import argparse
from Bio import SeqIO
from RouToolPa.Tools.Filter import Cookiecutter
from RouToolPa.Routines.Sequence import rev_com_generator
from RouToolPa.Routines.File import check_path, split_filename
from RouToolPa.Tools.Kmers import Jellyfish
from RouToolPa.Tools.Assemblers import MaSuRCA




parser = argparse.ArgumentParser()

parser.add_argument("-l", "--left_source_reads", action="store", dest="left_source_reads", required=True,
                    help="File with left source reads")
parser.add_argument("-r", "--right_source_reads", action="store", dest="right_source_reads", required=True,
                    help="File with right source reads")
parser.add_argument("-e", "--", action="store", dest="right_source_reads", required=True,
                    help="File with right source reads")
parser.add_argument("-v", "--right_source_reads", action="store", dest="right_source_reads", required=True,
                    help="File with right source reads")
parser.add_argument("-n", "--number_of_iterations", action="store", dest="number_of_iterations", type=int, default=20,
                    help="Number of iterations. Default: 20")
parser.add_argument("-i", "--initial_sequences", action="store", dest="initial_sequences",
                    help="Fasta with initial sequences")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-j", "--jellyfish_dir", action="store", dest="jellyfish_dir",
                    help="Jellyfish directory")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=int, default=1000000,
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L)")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers")
parser.add_argument("--mean_insert_size", action="store", dest="mean_insert_size", type=int, required=True,
                    help="Mean insert size")
parser.add_argument("--std_insert_size", action="store", dest="std_insert_size", type=int, required=True,
                    help="Standard deviation of insertsize")

parser.add_argument("-a", "--path_to_mavr", action="store", dest="path_to_mavr", default="./",
                    help="path_to_mavr")




parser.add_argument("-a", "--base_prefix", action="store", dest="base_prefix", default="jellyfish_db",
                    help="Name of kmer database. Default: jellyfish_db")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-b", "--count_both_strands", action="store_true", dest="count_both_strands",
                    help="Count kmers in both strands. NOTICE: only mer or its reverse-complement, whichever "
                         "comes first lexicographically, is stored and the count value is the number of "
                         "occurrences of both. So this option is not suitable for generating sets of forward "
                         "and reverse-complement kmers. For this case use -r/--add_reverse_complement option. "
                         "Not compatible with -r/--add_reverse_complement option.")
parser.add_argument("-r", "--add_reverse_complement", action="store_true", dest="add_rev_com",
                    help="Add reverse-complement sequences before counting kmers. "
                         "Works only for fasta sequences. "
                         "Not compatible with -b/--count_both_strands option")


parser = argparse.ArgumentParser()

args = parser.parse_args()
args.path_to_mavr = check_path(args.path_to_mavr)

MaSuRCA.threads = args.threads
Jellyfish.threads = args.threads
Jellyfish.path = args.jellyfish_path if args.jellyfish_path else ""

iteration_reference_file = args.initial_sequences
working_dir = os.getcwd()
abs_path_left_source_reads = os.path.abspath(args.left_source_reads)
abs_path_right_source_reads = os.path.abspath(args.right_source_reads)

"""
for filename in args.source_reads:
    ab
    if os.path.isabs(filename):
        abs_path_source_reads.append(filename)
    else:
        abs_path_source_reads.append("%s/%s" % (working_dir, filename))
"""

for iteration_index in range(1, args.number_of_iterations):

    os.chdir(working_dir)

    iteration = "iteration_%i" % iteration_index
    iteration_dir = "%s/%s" % (working_dir, iteration)
    iteration_ref = "%s/%s_reference.fasta" % (iteration_dir, iteration)
    iteration_ref_index = "%s/%s_reference.idx" % (iteration_dir, iteration)
    base_prefix = "%s/%s_reference_with_rev_com" % (iteration_dir, iteration)
    iteration_ref_with_rev_com = "%s/%s_reference_with_rev_com.fasta" % (iteration_dir, iteration)
    kmer_file = "%s_%i_mer.kmer" % (base_prefix, args.kmer_length)
    masurca_config_file = "masurca_%s.config" % iteration
    left_reads_prefix = split_filename(abs_path_left_source_reads)[1]
    right_reads_prefix = split_filename(abs_path_right_source_reads)[1]

    left_reads_se = "%s.se.fastq" % left_reads_prefix
    right_reads_se = "%s.se.fastq" % right_reads_prefix
    left_reads_filtered = "%s.filtered.fastq" % left_reads_prefix
    right_reads_filtered = "%s.filtered.fastq" % right_reads_prefix

    try:
        os.mkdir(iteration_dir)
    except OSError:
        pass

    shutil.copyfile(iteration_reference_file, iteration_ref)
    os.chdir(iteration_dir)
    iteration_reference_dict = SeqIO.index_db(iteration_ref_index, iteration_ref, format="fasta")
    SeqIO.write(rev_com_generator(iteration_reference_dict, yield_original_record=True),
                iteration_ref_with_rev_com, "fasta")

    Jellyfish.get_kmer_list(iteration_ref_with_rev_com, base_prefix, kmer_length=args.kmer_length,
                            hash_size=args.hash_size)
    Cookiecutter.extract(kmer_file, iteration_dir, abs_path_left_source_reads,
                         right_reads=abs_path_right_source_reads)
    os.system("%scripts/filter/restore_pairs.py -l %s -r %s -o %s" % (args.path_to_mavr,
                                                                      ",".join([left_reads_filtered, left_reads_se]),
                                                                      ",".join([right_reads_filtered, right_reads_se]),
                                                                      "iteration_%i_reads" % iteration_index))
    library = ("PE", args.mean_insert_size,
               args.std_insert_size,
               ["iteration_%i_reads_1.fastq" % iteration_index, "iteration_%i_reads_2.fastq" % iteration_index])

    MaSuRCA.generate_config([library], masurca_config_file, jellyfish_hash_size=args.hash_size, kmer_size="auto",
                            illumina_only_assembly=True, limit_jump_coverage=None, source="eukaryota",
                            trim_long_homopolymers=False, cgwErrorRate=0.15, ovlMemory="4GB",
                            minimum_count_kmer_in_error_correction=1)
    os.system(masurca_config_file)
    shutil.copyfile("../CA/10-gapclose/genome.scf.fasta", "iteration_%i_assembly.fasta" % iteration_index)



