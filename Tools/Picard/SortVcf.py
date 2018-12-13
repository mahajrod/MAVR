#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from Tools.Abstract import JavaTool


class SortVcf(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        #JavaTool.__init__(self, "CreateSequenceDictionary.jar", java_path=java_path, max_threads=max_threads,
        #                  jar_path=jar_path, max_memory=max_memory)
        JavaTool.__init__(self, "picard.jar SortVcf", java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=max_memory)

    def sort_vcf(self, input_file, output_file, seq_dict=None,
                 handling_mode="local",
                 max_memory_per_node=None,
                 job_name=None,
                 log_prefix=None,
                 error_log_prefix=None,
                 modules_list=None,
                 environment_variables_dict=None,
                 max_running_time=None):

        input_file_list = self.make_list_of_path_to_files_from_string(input_file) if isinstance(input_file, str) else self.make_list_of_path_to_files(input_file)

        options = ""

        for filename in input_file_list:
            options += " I=%s" % filename

        options += " O=%s" % output_file
        options += " SD=%s" % seq_dict if seq_dict else ""

        if handling_mode == 'local':
            self.execute(options=options)
        elif handling_mode == 'slurm':
            self.timelog = None
            slurm_cmd = self.execute(options=options, generate_cmd_string_only=True)

            job_id = self.slurm_run_job(job_name, log_prefix, slurm_cmd, error_log_prefix,
                                        "".join(self.split_filename(output_file)[:2]) + ".slurm",
                                        modules_list=modules_list,
                                        environment_variables_dict=environment_variables_dict,
                                        max_running_time=max_running_time,
                                        max_memory_per_node=max_memory_per_node)

            return job_id




if __name__ == "__main__":
    pass
