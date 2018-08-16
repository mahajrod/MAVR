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

    def parallel_align(self, database, query_fasta, output, split_dir="splited_input/",
                       splited_output_dir="splited_output_dir/",
                       database_type=None, query_type=None, add_header=True,
                       threads=None, remove_tmp_dirs=True,
                       async_run=False, external_process_pool=None):
        splited_dir = self.check_path(split_dir)
        splited_out_dir = self.check_path(splited_output_dir)
        self.safe_mkdir(splited_dir)
        self.safe_mkdir(splited_out_dir)

        common_options = self.parse_options(database_type=database_type,
                                            query_type=query_type,
                                            add_header=add_header)

        number_of_files = 5 * threads if threads else 5 * self.threads
        self.split_fasta(query_fasta, splited_dir, num_of_files=number_of_files)

        options_list = []
        input_list_of_files = sorted(os.listdir(splited_dir))
        output_file_list = []
        for filename in input_list_of_files:
            query = "%s/%s" % (splited_dir, filename)
            output_file = "%s/%s.out" % (splited_out_dir, filename)
            options = "%s" % database
            options += " %s" % query
            options += " %s" % common_options
            options += " %s" % output_file

            options_list.append(options)
            output_file_list.append(output_file)

        self.parallel_execute(options_list, cmd="blat", threads=threads, async_run=async_run,
                              external_process_pool=external_process_pool)

        cat_string = "cat %s > %s" % ("\t".join(output_file_list), output)
        os.system(cat_string)

        if remove_tmp_dirs:
            shutil.rmtree(splited_dir, ignore_errors=True)
