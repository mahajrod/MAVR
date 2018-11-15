#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from Tools.Abstract import Tool
from Parsers.Barrnap import CollectionBARRNAP


class BUSCO(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "run_BUSCO.py", path=path, max_threads=max_threads)
        self.kingdoms = ["euk", "bac", "mito", "arc"]

    def parse_options(self, input_fasta, label, busco_db, augustus_species, mode=None):

        options = " -c %i" % self.threads
        options += " -l %s" % busco_db
        options += " -m %s" % mode if mode else ""
        options += " -o %s" % label if label else ""
        options += " -sp %s" % augustus_species if augustus_species else ""
        options += " -i %s" % input_fasta

        return options

    def assess(self, input_fasta, label, busco_db, augustus_species, mode=None):

        options = self.parse_options(input_fasta, label, busco_db, augustus_species, mode=mode)

        self.execute(options=options)

    def assess_multiple_genomes(self, input_fasta_list, busco_db, augustus_species, label_list=None,
                                output_dir="./", mode=None):

        outdir = os.path.abspath(output_dir)

        for fasta, label in zip(input_fasta_list, label_list if label_list else ["A%i" % i for i in range(1,
                                                                                                          len(input_fasta_list) + 1)]):
            genome_fasta = os.path.abspath(fasta)
            options = self.parse_options(genome_fasta, label, busco_db, augustus_species, mode=mode)
            busco_dir = "%s/%s/" % (outdir, label)
            self.safe_mkdir(busco_dir)

            self.execute(options=options, directory=busco_dir)

if __name__ == "__main__":
    pass
