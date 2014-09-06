#!/usr/bin/env python
__author__ = 'mahajrod'

import os
from Tools.Picard import make_fasta_dict


def generate_mira_manifest_file(sample_name, reference_name, reference_file, read_group, fastq_file, out_file, technology, max_threads=5):
    with open(out_file, "w") as out_fd:
        out_fd.write(
        """
project = %s_mapping_assembly_to_%s
job = genome,mapping,accurate
parameters = -GE:not=%i -NW:mrnl=10000 -NW:cmrnl=warn -NW:cdrn=no -NW:cac=no -HS:nrr=50 -GE:kpmf=5 -SK:not=10
# The second part defines the sequencing data M# The data is logically divided into "readgroup# first, the reference sequence
readgroup
is_reference
data = %s
strain = %s
# now the Ion Torrent data
readgroup = %s
data = %s
technology = %s
strain = %s

""" % (sample_name, reference_name, max_threads, reference_file, reference_name, read_group, fastq_file, technology, sample_name))


workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq"
samples_list = ["001", "002", "003", "004"]

#reference_dir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/"
#reference_name = "LAN210_v0.10m"


reference_dir = "/home/mahajrod/genetics/desaminases/data/S288C_R64/"
reference_name = "S288C_R64"

reference = "%s%s.fasta" % (reference_dir, reference_name)
reference_dict = "%s%s.dict" % (reference_dir, reference_name)
reference_index = "%s%s" % (reference_dir, reference_name)

alignment_dir = "mira_alignment_%s" % reference_name
stat_file = "mira_stat_SNP_call_%s.t" % reference_name
intersect_dir = "mira_intersect_%s" % reference_name

os.chdir(reference_dir)
make_fasta_dict(reference, reference_dict, PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115/")
technology = "iontor"

for sample_name in samples_list:
    os.chdir(workdir)
    os.chdir(sample_name)
    os.system("mkdir -p %s" % alignment_dir)
    os.chdir(alignment_dir)
    mira_manifile = "%s.manifile" % sample_name
    for filename in os.listdir("../"):
        if filename[-6:] == ".fastq":
            fastq_file = filename
    fastq_file = "../" + fastq_file
    generate_mira_manifest_file(sample_name, reference_name, reference, "%s_reads" % sample_name, fastq_file, mira_manifile, technology)
    os.system("mira %s" % mira_manifile)
