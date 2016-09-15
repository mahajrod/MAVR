#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Sealer(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "abyss-sealer", path=path, max_threads=max_threads)

    def close_gaps(self, scaffold_fasta, forward_reads, reverse_reads, output_prefix,
                   kmer_size_list=(90, 80, 70, 60, 50, 40, 30), print_flanks=False, flank_len=None,
                   flank_distance=None, max_gap_len=None, bloom_size=None,
                   max_branches=None, dot_outfile=None, fix_errors=None, input_bloom_filter=None,
                   soft_mask_new_and_changed_bases=None, dont_soft_mask_bases=None, flank_mismatches=None,
                   max_mismatches=None, disable_all_limits=None, max_paths=None, trim_masked_from_ends=None,
                   do_not_discard_unchaste_reads=None, trim_quality=None):

        options = " -j %i" % self.threads
        options += " -S %s" %scaffold_fasta
        options += " -o %s" % output_prefix
        options += " %s" % (" ".join(map(lambda kmer_size: " -k%i" % kmer_size, kmer_size_list)))
        options += " --print-flanks" if print_flanks else ""
        options += " -L %i" % flank_len  if flank_len else ""
        options += " -D %i" % flank_distance if flank_distance else ""
        options += " -G %i" % max_gap_len if max_gap_len else ""
        options += " -b %s" % (bloom_size if isinstance(bloom_size, str) else "%i" % bloom_size) if bloom_size else ""
        options += " -B %i" % max_branches if max_branches else ""
        options += " -d %s" % dot_outfile if dot_outfile else ""
        options += " --fix-errors" if fix_errors else ""
        options += " -i %s" % input_bloom_filter if input_bloom_filter else ""
        options += " --mask" if soft_mask_new_and_changed_bases else ""
        options += " --no_mask" if dont_soft_mask_bases else ""
        options += " -m %i" % flank_mismatches if flank_mismatches else ""
        options += " -M %i" % max_mismatches if max_mismatches else ""
        options += " -P %i" % max_paths if max_paths else ""
        options += " -n" if disable_all_limits else ""
        options += " --trim-masked" if trim_masked_from_ends else ""
        options += " --no-chastity" if do_not_discard_unchaste_reads else ""
        options += " --trim-quality %i" % trim_quality if trim_quality else ""
        options += " %s" % forward_reads
        options += " %s" % reverse_reads

        self.execute(options)
        """
        Not implemented options

        --standard-quality       zero quality is `!' (33)
                               default for FASTQ and SAM files
        --illumina-quality       zero quality is `@' (64)
                               default for qseq and export files
        -r, --read-name=STR          only process reads with names that contain STR
        -s, --search-mem=N           mem limit for graph searches; multiply by the
                               number of threads (-j) to get the total mem used
                               for graph traversal [500M]
        -t, --trace-file=FILE        write graph search stats to FILE

        """

if __name__ == "__main__":
    pass
