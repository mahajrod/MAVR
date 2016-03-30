#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Quast(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "quast.py", path=path, max_threads=max_threads)

    def get_stat(self, contigs_file, out_dir="quast", reference=None, genes_file=None, operon_file=None, min_contig_length=None,
                 predict_genes=False, eukaryotic_genome=True, metagenome=False, scaffolds=False):
        options = " -T %i" % self.threads
        options += " -o %s" % out_dir
        options += " -R %s" % reference if reference else ""
        options += " -G %s" % genes_file if genes_file else ""
        options += " -O %s" % operon_file if operon_file else ""
        options += " --min-contig" % min_contig_length if min_contig_length else ""
        options += " -f" if predict_genes else ""
        options += " -e" if eukaryotic_genome else ""
        options += " -m" if metagenome else ""
        options += " -s" if scaffolds else ""

        options += " %s" % contigs_file
        self.execute(options)


if __name__ == "__main__":
    pass