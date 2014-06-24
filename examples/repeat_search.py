#!/usr/bin/env python
import os
from Bio import SeqIO

#sys.path.append("../")
from Pipelines.RepeatSearch import *


seq_dict = SeqIO.to_dict(SeqIO.parse("data/Acipenseriformes_complete_MH_genome.gb", "genbank"))
workdir = os.getcwd()

for seq_id in seq_dict:
    os.chdir(workdir)
    os.system("mkdir -p %s" % seq_id)
    os.chdir(seq_id)
    SeqIO.write([seq_dict[seq_id]], seq_id + ".fasta", "fasta")
    make_blast_db(seq_id + ".fasta", seq_id)
    sliding_window_blastn(seq_dict[seq_id].seq, seq_id)