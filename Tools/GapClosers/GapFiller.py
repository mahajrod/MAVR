#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class GapFiller(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "SSPACE_Standard_v3.0.pl", path=path, max_threads=max_threads)

    def close_gaps(self, library_file, scaffold_fasta, output_basename, min_overlap_len=None,
                   min_read_number_for_base_calling=None, read_with_extention_fraction_for_closing=None,
                   max_gap_len_difference_to_stop=None, num_reads_to_trim=None, iteration_number=None,
                   min_contig_overlap=None, max_gaps_for_bowtie=None, skip_read_processing=None):

        options = " -T %i" % self.threads
        options += " -l %s" % library_file
        options += " -s %s" % scaffold_fasta

        options += " -m %i" % min_overlap_len if min_overlap_len else ""
        options += " -o %i" % min_read_number_for_base_calling if min_read_number_for_base_calling else ""
        options += " -r %f" % read_with_extention_fraction_for_closing if read_with_extention_fraction_for_closing else ""
        options += " -d %i" % max_gap_len_difference_to_stop if max_gap_len_difference_to_stop else ""
        options += " -t %i" % num_reads_to_trim if num_reads_to_trim else ""
        options += " -n %i" % min_contig_overlap if min_contig_overlap else ""
        options += " -i %i" % iteration_number if iteration_number else ""
        options += " -g %i" % max_gaps_for_bowtie if max_gaps_for_bowtie else ""
        options += " -S" if skip_read_processing else ""
        options += " -b %s" % output_basename if output_basename else ""

        self.execute(options)

if __name__ == "__main__":
    pass


