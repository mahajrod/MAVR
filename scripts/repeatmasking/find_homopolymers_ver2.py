#!/usr/bin/env python3
__author__ = 'mahajrod'
# !!!!!!!!!! Input is normal fasta from flybase
import argparse
import os
import numpy as np
import matplotlib

matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()
from Bio import SeqIO

from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="in_fasta",
                    help="input fasta file with sequences in one line")
parser.add_argument("-d", "--index_file", action="store", dest="idx_fasta", default="temp.idx",
                    help="Biopython index of input fasta (will be created is abbsent)")
parser.add_argument("-o", "--output", action="store", dest="out_file",
                    help="output file")
parser.add_argument("-n", "--nucleotide", action="store", dest="nucleotide",
                    help="nucleotide to search homopolymers")
parser.add_argument("-m", "--min_size", action="store", dest="min_size", type=int, default=5,
                    help="minimum size of homopolymer to be counted in")
parser.add_argument("-s", "--search_type", action="store", dest="search_type", default="perfect",
                    help="type of search, possible values: 'perfect' or 'non_perfect'")
parser.add_argument("-x", "--max_single_insert_size", action="store",
                    dest="max_single_insert_size", type=int, default=1,
                    help="maximum allowed size of single insertion. is not taken into account is 'perfect' search type is selected")
parser.add_argument("-y", "--max_number_of_insertions", action="store",
                    dest="max_number_of_insertions", type=int, default=2,
                    help="maximum number of single insertions. is not taken into account is 'perfect' search type is selected")
parser.add_argument("-z", "--max_total_insert_length", action="store",
                    dest="max_total_insert_length", type=int, default=None,
                    help="maximum total size of insertions. is not taken into account is 'perfect' search type is selected")
parser.add_argument("-p", "--hist_prefix", action="store", dest="hist_prefix", default="histogram",
                    help="prerefix of histogram images")
parser.add_argument("-e", "--species", action="store", dest="species", default="Unknown",
                    help="species")
parser.add_argument("-r", "--seq_type", action="store", dest="seq_type", default="Unknown_seq",
                    help="Type of sequences. string, used in histogram title")
parser.add_argument("-l", "--id_file", action="store", dest="id_file", default="id_file_list.t",
                    help="File with ids of UTRS with homopolymers")
parser.add_argument("-g", "--gene_id_file", action="store", dest="gene_id_file", default="gene_id_file_list.txt",
                    help="File with ids of genes with homopolymers in UTRs")
parser.add_argument("-a", "--gene_name_file", action="store", dest="gene_name_file", default="gene_name_file_list.txt",
                    help="File with names of genes with homopolymers in UTRs")
parser.add_argument("-f", "--outfiltered_gene_name_file", action="store", dest="outfiltered_gene_name_file", default="outfiltered_gene_name_file_list.txt",
                    help="File with names of outfiltered genes")
args = parser.parse_args()

UTRs_dict = SeqIO.index_db(args.idx_fasta, args.in_fasta, format="fasta")

max_length_list = []
number_poly_list = []
number_of_UTRs = 0
id_fd = open(args.id_file, "w")
gene_id_fd = open(args.gene_id_file, "w")
gene_name_fd = open(args.gene_name_file, "w")
out_filtered_genes = []
out_filterd_gene_name_fd = open(args.outfiltered_gene_name_file, "w")
with open(args.out_file, "w") as out_fd:
    out_fd.write("#main_id\tUTR_length\tmax_homopolymer_length\tnumber_homopolymers\tCoordinates_list\tOther_ids\n")
    for record_id in UTRs_dict:
        sequence = UTRs_dict[record_id].seq
        number_of_UTRs += 1

        coords_list, length_list = SequenceRoutines.find_homopolymers(sequence, args.nucleotide, min_size=args.min_size,
                                                     search_type=args.search_type,
                                                     max_single_insert_size=args.max_single_insert_size,
                                                     max_total_insert_length=args.max_total_insert_length,
                                                     max_number_of_insertions=args.max_number_of_insertions)
        gene_name = UTRs_dict[record_id].description.split(";")[2].split("=")[1].split("-")[0]
        if not coords_list:
            out_filterd_gene_name_fd.write(gene_name + "\n")
            continue
        id_fd.write(record_id + "\n")
        gene_id = UTRs_dict[record_id].description.split(";")[5].split("=")[1]   # extract parent id from description of flybase fasta
        gene_id_fd.write(gene_id + "\n")
          # extract gene name from description of flybase fasta
        gene_name_fd.write(gene_name + "\n")
        max_length = max(length_list)
        number_of_homopolymers = len(length_list)
        max_length_list.append(max_length)
        number_poly_list.append(number_of_homopolymers)
        coords_str_list = map(lambda x: "(%i,%i)" % (x[0], x[1]), coords_list)
        out_fd.write("%s\t%i\t%i\t%i\t%s\n" % (record_id, len(sequence), max_length, number_of_homopolymers,
                                               ",".join(coords_str_list)))

for file_fd in id_fd, gene_id_fd, gene_name_fd, out_filterd_gene_name_fd:
    file_fd.close()

max_homopolymer = max(max_length_list)
bins = np.linspace(args.min_size, max_homopolymer, max_homopolymer - args.min_size + 1)

plt.figure(1, figsize=(6, 6))
plt.subplot(111)
plt.hist(max_length_list, bins=bins, label="Totaly %s: %i\nWith poly(T) %i+: %i\nLongest poly(%s): %i bp" %
                                           (args.seq_type,
                                            number_of_UTRs,
                                            args.min_size,
                                            len(max_length_list),
                                            args.nucleotide,
                                            max_homopolymer))
plt.xlabel("Size of longest poly(%s)" % args.nucleotide)
plt.ylabel("Number of %s" % args.seq_type)
plt.title("Poly(%s) in %s of %s" % (args.nucleotide, args.seq_type, args.species))
plt.legend()
for extension in ".svg", ".png", ".eps":
    plt.savefig(args.hist_prefix + extension)
plt.close()

max_number_of_polyT = max(number_poly_list)
xbins = bins
ybins = bins = np.linspace(1, max_number_of_polyT, max_number_of_polyT)
plt.figure(2, figsize=(6, 6))
plt.subplot(111)
plt.hist2d(max_length_list, number_poly_list, bins=(xbins, ybins), cmin=1)
plt.xlim(xmax=max_homopolymer)
plt.ylim(ymax=max_number_of_polyT)
plt.xlabel("Size of longest poly(%s)" % args.nucleotide)
plt.ylabel("Number of poly(%s)" % args.nucleotide)
plt.title("Poly(%s) in %s of %s" % (args.nucleotide, args.seq_type, args.species))
plt.legend()
for extension in ".svg", ".png", ".eps":
    plt.savefig(args.hist_prefix + "_2d" + extension)
plt.close()

os.remove(args.idx_fasta)
