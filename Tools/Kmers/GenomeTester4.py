#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

from Tools.Abstract import Tool


class Glistmaker(Tool):
    """ """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "glistmaker", path=path, max_threads=max_threads)

    def parse_options(self, in_file_list, output_file, kmer_length=16, frequency_cutoff=None, threads=None,
                      max_tmp_table_number=None, max_tmp_table_size=None):

        options = " %s" % (in_file_list if isinstance(in_file_list, str) else " ".join(in_file_list))
        options += " -o %s" % output_file if output_file else ""
        options += " -w %i" % kmer_length
        options += " -c %i" % frequency_cutoff if frequency_cutoff else ""
        options += " --num_threads %i" % (threads if threads else self.threads)
        options += " --max_tables %i" % max_tmp_table_number if max_tmp_table_number else ""
        options += " --table_size %i" % max_tmp_table_size if max_tmp_table_size else ""

        return options

    def count(self, in_file_list, output_file, kmer_length=16, frequency_cutoff=None, threads=None,
              max_tmp_table_number=None, max_tmp_table_size=None):

        options = self.parse_options(in_file_list, output_file, kmer_length=kmer_length,
                                     frequency_cutoff=frequency_cutoff, threads=threads,
                                     max_tmp_table_number=max_tmp_table_number,
                                     max_tmp_table_size=max_tmp_table_size)

        self.execute(options)

    def generate_kmer_lists_for_primer3(self, in_file_list, species_prefix, threads=None,
                                        max_tmp_table_number=None, max_tmp_table_size=None):
        #out_kmer_16 = "%s_16.list" % species_prefix
        #out_kmer_11 = "%s_11.list" % species_prefix

        options_kmer_16 = self.parse_options(in_file_list, species_prefix, kmer_length=16, threads=threads,
                                             max_tmp_table_number=max_tmp_table_number,
                                             max_tmp_table_size=max_tmp_table_size)
        options_kmer_11 = self.parse_options(in_file_list, species_prefix, kmer_length=11, threads=threads,
                                             max_tmp_table_number=max_tmp_table_number,
                                             max_tmp_table_size=max_tmp_table_size)

        self.execute(options_kmer_16)
        self.execute(options_kmer_11)
