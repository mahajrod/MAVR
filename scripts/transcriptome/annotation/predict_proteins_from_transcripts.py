#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.BLAST import BLASTp
from RouToolPa.Tools.HMMER import HMMER3
from RouToolPa.Tools.Annotation import TransDecoder
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Fasta file with transcripts")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output fasta file with complete peptides")
parser.add_argument("-g", "--genetic_code", action="store", dest="genetic_code",
                    help="Genetic code to use. Possible genetic codes:"
                         "\tuniversal (default)\n"
                         "\tEuplotes\n"
                         "\tTetrahymena\n"
                         "\tCandida\n"
                         "\tAcetabularia\n"
                         "\tMitochondrial-Canonical\n"
                         "\tMitochondrial-Vertebrates\n"
                         "\tMitochondrial-Arthropods\n"
                         "\tMitochondrial-Echinoderms\n"
                         "\tMitochondrial-Molluscs\n"
                         "\tMitochondrial-Ascidians\n"
                         "\tMitochondrial-Nematodes\n"
                         "\tMitochondrial-Platyhelminths\n"
                         "\tMitochondrial-Yeasts\n"
                         "\tMitochondrial-Euascomycetes\n"
                         "\tMitochondrial-Protozoans")
parser.add_argument("-m", "--min_protein_length", action="store", dest="min_prot_len", default=100, type=int,
                    help="Minimum protein length. Default: 100")
parser.add_argument("-S", "--analyze_only_top_strand", action="store_true", dest="analyze_only_top_strand",
                    help="Analyze only top strand. Default: False")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Minimum protein length. Default: 100")
parser.add_argument("-p", "--pfam_database", action="store", dest="pfam_database",
                    help="File with hmm3 profiles from pfam database")
parser.add_argument("-b", "--blast_database", action="store", dest="blast_database",
                    help="Blast database with domains")
parser.add_argument("-r", "--min_orf_len_if_no_other_evidence", action="store", type=int,
                    dest="min_orf_len_if_no_other_evidence",
                    help="Retain all ORFs found that are equal or longer than these many nucleotides "
                         "even if no other evidence")
parser.add_argument("-a", "--file_with_orfs_for_training", action="store", dest="file_with_orfs_for_training",
                    help="FASTA file with ORFs to train Markov Mod for protein identification; "
                         "otherwise longest non-redundant ORFs used")
parser.add_argument("-n", "--number_of_top_orfs_for_training", action="store", type=int,
                    dest="number_of_top_orfs_for_training",
                    help="If no --number_of_top_orfs_for_training, top longest ORFs to train Markov "
                         "Model (hexamer stats) (default: 500)")
parser.add_argument("-c", "--hmmer_dir", action="store", dest="hmmer_dir", default="",
                    help="Directory with hmmer v3.1 binaries")
parser.add_argument("-d", "--blast_dir", action="store", dest="blast_dir", default="",
                    help="Directory with BLAST+ binaries")

args = parser.parse_args()

input_filename_list = FileRoutines.split_filename(args.input)
input_filename = input_filename_list[1] + input_filename_list[2]

workdir_dir = "%s.transdecoder_dir/" % input_filename
pep_from_longest_orfs = "%s/longest_orfs.pep" % workdir_dir

hmmscan_dir = "hmmscan_vs_pfam/"
blastp_dir = "blastp_vs_uniref/"

FileRoutines.safe_mkdir(hmmscan_dir)
FileRoutines.safe_mkdir(blastp_dir)

hmmscan_splited_fasta_dir = "%ssplited_fasta_dir/" % hmmscan_dir
splited_domtblout_dir = "%ssplited_domtblout_dir/" % hmmscan_dir

hmmscan_vs_pfam_prefix = "%s.pfam" % input_filename
hmmscan_vs_pfam_output = "%s/%s.hits" % (hmmscan_dir, hmmscan_vs_pfam_prefix)
domtblout_outfile = "%s/%s.domtblout" % (hmmscan_dir, hmmscan_vs_pfam_prefix)

blastp_outfile = "%s%s.blastp.hits" % (blastp_dir, input_filename) if args.blast_database else None
blastp_split_dir = "%ssplited_fasta_dir/" % blastp_dir
blastp_splited_output_dir = "%ssplited_output_dir" % blastp_dir
HMMER3.path = args.hmmer_dir
HMMER3.threads = args.threads
BLASTp.path = args.blast_dir
BLASTp.threads = args.threads

TransDecoder.extract_longest_orfs(args.input, genetic_code=args.genetic_code,
                                  analyze_only_top_strand=args.analyze_only_top_strand,
                                  minimum_protein_length=args.min_prot_len)
if args.pfam_database:
    HMMER3.parallel_hmmscan(args.pfam_database, pep_from_longest_orfs, hmmscan_vs_pfam_prefix,
                            hmmscan_dir, dont_output_alignments=True)
if args.blast_database:
    BLASTp.parallel_blastp(pep_from_longest_orfs, args.blast_database, outfile=blastp_outfile,
                           evalue=0.00001, output_format=6, blast_options=" -max_target_seqs 1",
                           combine_output_to_single_file=True, split_dir=blastp_split_dir,
                           splited_output_dir=blastp_splited_output_dir)

TransDecoder.predict_pep(args.input, pfam_hits=domtblout_outfile,
                         blastp_hits=blastp_outfile,
                         minimum_orf_length_if_no_other_evidence=args.min_orf_len_if_no_other_evidence,
                         file_with_orfs_for_training=args.file_with_orfs_for_training,
                         number_of_top_orfs_for_training=args.number_of_top_orfs_for_training)
