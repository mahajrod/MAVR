#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from Tools.Abstract import Tool
from Parsers.Barrnap import CollectionBARRNAP

from Routines.File import split_filename, save_mkdir
from CustomCollections.GeneralCollections import SynDict, IdList


class Barrnap(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "barrnap", path=path, max_threads=max_threads)
        self.kingdoms = ["euk", "bac", "mito", "arc"]

    def parse_options(self, input_fasta, kingdom, gff_file, log_file, length_cutoff=None, reject_cutoff=None,
                      evalue_cutoff=None):

        if kingdom not in self.kingdoms:
            raise ValueError("Wrong kingdom of life. Allowed: euk, bac, mito, arc")

        options = " --threads %i" % self.threads
        options += " --kingdom" % kingdom
        options += " --lencutoff %f" if length_cutoff else ""
        options += " --reject %f" if reject_cutoff else ""
        options += " --evalue %f" if evalue_cutoff else ""
        options += " %s" % input_fasta
        options += " > %s " % gff_file
        options += " 2> %s" % log_file

        return options

    def annotate_rrna(self, input_fasta, kingdom, output_prefix, length_cutoff=None, reject_cutoff=None,
                      evalue_cutoff=None):

        gff_file = "%s.%s.gff" % (output_prefix, kingdom)
        log_file = "%s.%s.log" % (output_prefix, kingdom)
        options = self.parse_options(input_fasta, kingdom, gff_file, log_file, length_cutoff=length_cutoff,
                                     reject_cutoff=reject_cutoff, evalue_cutoff=evalue_cutoff)
        self.execute(options)
        self.analyze_report(log_file, "%s.%s" % (output_prefix, kingdom))

        if kingdom == "euk":
            mito_gff_file = "%s.mito.gff" % output_prefix
            mito_log_file = "%s.mito.log" % output_prefix
            mito_options = self.parse_options(input_fasta, "mito", mito_gff_file, mito_log_file,
                                              length_cutoff=length_cutoff,
                                              reject_cutoff=reject_cutoff, evalue_cutoff=evalue_cutoff)

            self.execute(mito_options)
            self.analyze_report(mito_log_file, "%s.mito" % output_prefix)

    @staticmethod
    def analyze_report(report_file, output_prefix):
        count_file = "%s.chrom.count" % output_prefix
        total_file = "%s.total.count" % output_prefix

        report = CollectionBARRNAP(from_file=True, barrnap_file=report_file)

        count_dict, total_count_dict = report.count_types(output_file=count_file,
                                                          total_output_file=total_file,
                                                          return_mode="both")
        return count_dict, total_count_dict


if __name__ == "__main__":
    pass
