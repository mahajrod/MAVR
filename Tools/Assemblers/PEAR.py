#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class PEAR(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "pear", path=path, max_threads=max_threads)

    def get_stat(self, forward_reads, reverse_reads, output, p_value=None, min_overlap=None, max_assembly_length=None,
                 min_assembly_length=None, min_trim_length=None, quality_threshold=None, max_uncalled_base=None,
                 test_method=None, disable_empirical_base_frequencies=False, score_method=None, phred_base=None,
                 memory="4G", max_quality_for_resulting_scores=0, threads=1, replace_non_equal_overlapping_bases_by_N=None):

        options = " --forward-fastq %s" % forward_reads
        options += " --reverse-fastq %s" % reverse_reads
        options += " --output %s" % output
        options += " --p-value %f" % p_value if p_value else ""
        options += " --min-overlap %i" % min_overlap if min_overlap else ""
        options += " --max-assembly-length %i" % max_assembly_length if max_assembly_length else ""
        options += " --min-assembly-length %i" % min_assembly_length if min_assembly_length else ""
        options += " --min-trim-length %i" % min_trim_length if min_trim_length else ""

        options += " --quality-threshold %i" % quality_threshold if quality_threshold else ""
        options += " --max-uncalled-base %f" % max_uncalled_base if max_uncalled_base else ""
        options += " --test-method %i" % test_method if test_method else ""
        options += " --empirical-freqs" if disable_empirical_base_frequencies else ""
        options += " --score-method %i" % score_method if score_method else ""
        options += " --phred-base %i" % phred_base if phred_base else ""
        options += " --memory %s" % memory if memory else ""
        options += " --cap %i" % max_quality_for_resulting_scores if max_quality_for_resulting_scores else ""

        options += " --threads %i" % threads if threads else ""
        options += " --nbase" if replace_non_equal_overlapping_bases_by_N else ""
        self.execute(options)


if __name__ == "__main__":
    pass
