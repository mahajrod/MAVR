#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
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
    '''listens for messages on the q, writes to file. '''
    #print (queue, prot_fd, gene_fd)
    protein_fd = open(args.prefix + "_proteins.fam", "w")
    genes_fd = open(args.prefix + "_genes.fam", "w")
    while 1:
        m = queue.get()
        #print m
        if m == 'kill':
            break
        #print m
        protein_fd.write("%s\t%s\n" % (m[0], ",".join(m[2])))
        genes_fd.write("%s\t%s\n" % (m[0], ",".join(m[1])))
        protein_fd.flush()
        genes_fd.flush()
    protein_fd.close()
    genes_fd.close()

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



