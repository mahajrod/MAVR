#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Trimmer(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "V2_trim.exe", path=path, max_threads=max_threads)

    def parse_options(self, reads_prefix, kmer_file, output_prefix, write_fasta=False, input_format="fastq"):

        options = " %s" % reads_prefix
        options += " %s" % output_prefix
        options += " %i" % self.threads
        options += " %i" % (1 if write_fasta else 0)
        options += " %s" % input_format
        options += " %s" % kmer_file
        #options += " > %s.stats 2>&1" % output_prefix

        return options

    def trim_adapters(self, reads_prefix, adapter_kmer_file, output_prefix, write_fasta=False, input_format="fastq"):

        options = self.parse_options(reads_prefix, adapter_kmer_file, output_prefix,
                                     write_fasta=write_fasta, input_format=input_format)

        self.execute(options=options)

    def mass_trim(self, sample_dir, output_dir, adapter_kmer_file, sample_list=None,
                  write_fasta=False, input_format="fastq"):

        samples = self.get_sample_list(sample_dir, sample_list)

        self.safe_mkdir(output_dir)

        for sample in samples:
            sample_read_prefix = "%s/%s/%s" % (sample_dir, sample, sample)
            out_sample_dir = "%s/%s/" % (output_dir, sample)
            out_sample_prefix = "%s/%s" % (out_sample_dir, sample)
            self.safe_mkdir(out_sample_dir)

            self.trim_adapters(sample_read_prefix, adapter_kmer_file, out_sample_prefix,
                               write_fasta=write_fasta, input_format=input_format)


if __name__ == "__main__":
    pass


