#!/usr/bin/env python
import os
import shutil
from Tools.Abstract import Tool


class BLAT(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blat", path=path, max_threads=max_threads)

    def parse_options(self, database_type=None, query_type=None, add_header=True):

        options = ""
        options += " -t=%s" % database_type if database_type else ""
        options += " -q=%s" % query_type if query_type else ""
        options += " -noHead" if not add_header else ""

        return options

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
        if add_header:
            with open(output, "w") as out_fd:
                pass
            os.system("cat %s > %s" % (output_file_list[0], output))
            for filename in output_file_list[1:]:
                sed_string = "sed -n 6,$p %s > %s" % (filename, output)
                os.system(sed_string)
        else:
            cat_string = "cat %s > %s" % ("\t".join(output_file_list), output)
            os.system(cat_string)

        if remove_tmp_dirs:
            shutil.rmtree(splited_dir, ignore_errors=True)
            shutil.rmtree(splited_output_dir, ignore_errors=True)

    def align_pep_for_scaffolding(self, database, query_fasta, output, split_dir="splited_input/",
                                  splited_output_dir="splited_output_dir/",
                                  threads=None, remove_tmp_dirs=True,
                                  async_run=False, external_process_pool=None):

        self.parallel_align(database, query_fasta, output, split_dir=split_dir,
                            splited_output_dir=splited_output_dir,
                            database_type="dnax", query_type="prot", add_header=False,
                            threads=threads, remove_tmp_dirs=remove_tmp_dirs,
                            async_run=async_run, external_process_pool=external_process_pool)

    def align_transcripts_for_scaffolding(self, database, query_fasta, output, split_dir="splited_input/",
                                          splited_output_dir="splited_output_dir/",
                                          threads=None, remove_tmp_dirs=True,
                                          async_run=False, external_process_pool=None):

        self.parallel_align(database, query_fasta, output, split_dir=split_dir,
                            splited_output_dir=splited_output_dir,
                            database_type="dna", query_type="dna", add_header=False,
                            threads=threads, remove_tmp_dirs=remove_tmp_dirs,
                            async_run=async_run, external_process_pool=external_process_pool)