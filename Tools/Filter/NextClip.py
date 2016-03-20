#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class NextClip(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "nextclip", path=path, max_threads=max_threads)

    def filter(self, forward_reads, reverse_reads, output_prefix, log_file, min_len=25, number_of_reads=None,
               trim_ends=0, remove_duplicates=True, use_category_e=False, duplicates_log=False, strict_match=None,
               relaxed_match=None,  adaptor_sequence=None):
        # strict_match and relaxed_match have to be tuples with 2 numbers
        # Default values:
        # [-x | --strict_match] Strict alignment matches (default '34,18')
        # [-y | --relaxed_match] Relaxed alignment matches (default '32,17')

        options = " --input_one %s" % forward_reads
        options += " --input_two %s" % reverse_reads
        options += " --output_prefix %s" % output_prefix
        options += " --log %s" % log_file
        options += " --min_length %i" % min_len
        options += " --number_of_reads %i" % number_of_reads if number_of_reads else ""
        options += " --trim_ends %i" % trim_ends
        options += " --remove_duplicates" if remove_duplicates else ""
        options += " --use_category_e" if use_category_e else ""
        options += " --duplicates_log %s" % duplicates_log if duplicates_log else ""
        options += " --strict_match %i,%i" % (strict_match[0], strict_match[1]) if strict_match else ""
        options += " --relaxed_match %i,%i" % (relaxed_match[0], relaxed_match[1]) if relaxed_match else ""
        options += " --adaptor_sequence %s" % adaptor_sequence if adaptor_sequence else ""

        self.execute(options, cmd="filter_by_mean_quality")


if __name__ == "__main__":
    pass



