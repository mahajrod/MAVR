#!/usr/bin/env python
import os
import shutil
from Tools.Abstract import Tool


class GEM(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "gem", path=path, max_threads=max_threads)

    def parse_options(self, reference=None, threads=None, data_type=None, output_prefix=None,
                      kmer_length=None, gem_index=None):

        options = " -T %i" % (threads if threads else self.threads)
        options += " -c %s" % data_type if data_type else ""
        options += " -i %s" % reference if reference else ""
        options += " -o %s" % output_prefix if output_prefix else ""
        options += " -l %i" % kmer_length if kmer_length else ""
        options += " -I %s" % gem_index if gem_index else ""

        return options

    def index(self, reference, output_prefix, threads=None, data_type="dna"):

        options = self.parse_options(reference=reference,
                                     threads=threads,
                                     data_type=data_type,
                                     output_prefix=output_prefix)

        self.execute(options=options, cmd="gem-indexer")

    def create_map_of_mappability(self, kmer_list, output_prefix, gem_index=None, reference=None, threads=None,
                                  convert_to_wig=True):
        if (gem_index is None) and (reference is None):
            raise ValueError("ERROR!!! Neither Gem index nor reference was set!")

        if gem_index is None:
            self.index(reference, output_prefix, threads, data_type='dna')

        gem_idx = gem_index if gem_index else "%s.gem" % output_prefix

        for kmer_len in [kmer_list] if isinstance(kmer_list, int) else kmer_list:
            kmer_mappability_map_prefix = "%s.k%i" % (output_prefix, kmer_len)
            options = self.parse_options(threads=threads, output_prefix=kmer_mappability_map_prefix,
                                         gem_index=gem_idx, kmer_length=kmer_len)

            self.execute(options=options, cmd="gem-mappability")
            if convert_to_wig:
                self.convert_mappability_map_to_wig(map_of_mappability="%s.mappability" % kmer_mappability_map_prefix,
                                                    output_prefix=kmer_mappability_map_prefix,
                                                    gem_index=gem_idx)

    def convert_mappability_map_to_wig(self, map_of_mappability, output_prefix, threads=None, gem_index=None, reference=None):

        if (gem_index is None) and (reference is None):
            raise ValueError("ERROR!!! Neither Gem index nor reference was set!")

        if gem_index is None:
            self.index(reference, output_prefix, threads, data_type='dna')

        gem_idx = gem_index if gem_index else "%s.gem" % output_prefix
        options = self.parse_options(gem_index=gem_idx, output_prefix=output_prefix)
        options += " -i %s" % map_of_mappability

        self.execute(options=options, cmd="gem-2-wig")

