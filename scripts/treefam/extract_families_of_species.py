#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
Tested on treefam database 9.
Species with problems during extracting:
    UNSOLVED: uncommon format in treefam database(even ids were not extracted), solution - write custom scripts for extraction
        bursaphelenchus_xylophilus
        capitella_teleta
        helobdella_robusta
        heterorhabditis_bacteriophora
        lottia_gigantea
        meloidogyne_hapla
        monosiga_brevicollis
        proterospongia_sp
        strongyloides_ratti
    SOLVED: different ids in sequence and families (presence ":" in ids in files with families - in sequence files is replaced by "_"), solved by replacement in files with extracted ids
        tribolium_castaneum
        amphimedon_queenslandica
        schizosaccharomyces_pombe
"""

import os
import sys
import argparse

import multiprocessing as mp
import subprocess as sb


def extract_species_genes(family_file_name, queue):
    family_name = family_file_name[:-8]
    gene_lines = sb.Popen("grep '%s' %s/%s | tee %s/%s"
                          % (args.species, args.input, family_file_name, args.output, family_file_name),
                          shell=True, stdout=sb.PIPE)
    protein_list = []
    gene_list = []
    for line in gene_lines.stdout:
        line_list = line.strip().split()
        protein_list.append(line_list[2])
        gene_list.append(line_list[7])
    if gene_list:
        queue.put((family_name, gene_list, protein_list))
        return family_name, gene_list, protein_list


def listener(queue):
    """listens for messages on the queue, writes to file."""
    protein_fd = open(args.prefix + "_proteins.fam", "w")
    genes_fd = open(args.prefix + "_genes.fam", "w")
    protein_ids_fd = open(args.prefix + "_in_treefam_families_protein.ids", "w")
    gene_ids_fd = open(args.prefix + "_in_treefam_families_genes.ids", "w")
    while 1:
        m = queue.get()
        if m == 'kill':
            break
        protein_fd.write("%s\t%s\n" % (m[0], ",".join(m[2])))
        genes_fd.write("%s\t%s\n" % (m[0], ",".join(m[1])))
        for protein_id in m[2]:
            protein_ids_fd.write(protein_id + "\n")
        for gene_id in m[1]:
            gene_ids_fd.write(gene_id + "\n")
        protein_fd.flush()
        genes_fd.flush()
        protein_ids_fd.flush()
        gene_ids_fd.flush()

    protein_fd.close()
    genes_fd.close()
    protein_ids_fd.close()
    gene_ids_fd.close()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_directory", action="store", dest="input",
                    help="Directory with TreeFam database.")
parser.add_argument("-s", "--species", action="store", dest="species",
                    help="Latin name of species to extract. Example: mus_musculus, homo_sapiens etc")
parser.add_argument("-o", "--output_directory", action="store", dest="output",
                    help="Output directory. Default: species")
parser.add_argument("-f", "--prefix", action="store", dest="prefix",
                    help="Prefix of files with families.")
parser.add_argument("-p", "--threads", action="store", dest="threads", type=int,
                    help="Number of threads. Default: 1")
args = parser.parse_args()

if args.output is None:
    args.output = args.species

if args.prefix is None:
    args.prefix = args.species
try:
    os.mkdir(args.output)
except OSError:
    pass

list_of_files = os.listdir(args.input)

nhx_files = []
for entry in list_of_files:
    #print entry
    if entry[-8:] == '.nhx.emf':
        nhx_files.append(entry)
print("Found %i families" % len(nhx_files))
print("Starting extraction...")

manager = mp.Manager()
queue = manager.Queue()
process_pool = mp.Pool(args.threads)

watcher = process_pool.apply_async(listener, (queue, ))

jobs = []

for entry in nhx_files:
    #print entry
    job = process_pool.apply_async(extract_species_genes, (entry, queue))
    jobs.append(job)

for job in jobs:
    #print job
    job.get()

    #now we are done, kill the listener
queue.put('kill')
process_pool.close()



