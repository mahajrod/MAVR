#!/usr/bin/env python
import os
from Tools.Abstract import Tool

from Routines import FileRoutines

class PRANK(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "prank", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(tree_file=None, output_format=None, show_xml=None, show_tree=None,
                             show_ancestral_sequences=None, show_evolutionary_events=None,
                             showall=None, compute_posterior_support=None, njtree=None, skip_insertions=False,
                             codon_alignment=None, translated_alignment=None):
        # TODO: add rest of options
        options = " -t=%s" % tree_file if tree_file else ""
        options += " -f=%s" % output_format if output_format else ""
        options += " -showxml" if show_xml else ""
        options += " -showtree" if show_tree else ""
        options += " -showanc" if show_ancestral_sequences else ""
        options += " -showevents" if show_evolutionary_events else ""
        options += " -showall" if showall else ""
        options += " -support" if compute_posterior_support else ""
        options += " -njtree" if njtree else ""
        options += " -F" if skip_insertions else ""
        options += " -codon" if codon_alignment else ""
        options += " -translate" if translated_alignment else ""

        return options

    def align(self, sequence_file, output, tree_file=None, output_format=None, show_xml=None,
              show_tree=None, show_ancestral_sequences=None, show_evolutionary_events=None,
              showall=None, compute_posterior_support=None, njtree=None, skip_insertions=False,
              codon_alignment=None, translated_alignment=None):
        # TODO: add rest of options
        options = " -d=%s" % sequence_file
        options += " -o=%s" % output

        options += self.parse_common_options(tree_file=tree_file, output_format=output_format, show_xml=show_xml,
                                             show_tree=show_tree, show_ancestral_sequences=show_ancestral_sequences,
                                             show_evolutionary_events=show_evolutionary_events, showall=showall,
                                             compute_posterior_support=compute_posterior_support, njtree=njtree,
                                             skip_insertions=skip_insertions, codon_alignment=codon_alignment,
                                             translated_alignment=translated_alignment)
        self.execute(options)

    def parallel_align(self, list_of_files, output_directory, output_suffix=None, tree_file=None, output_format=None, show_xml=None,
                       show_tree=None, show_ancestral_sequences=None, show_evolutionary_events=None,
                       showall=None, compute_posterior_support=None, njtree=None, skip_insertions=False,
                       codon_alignment=None, translated_alignment=None,
                       cmd_log_file=None,
                       cpus_per_task=1,
                       handling_mode="local",
                       job_name=None,
                       log_prefix=None,
                       error_log_prefix=None,
                       max_jobs=None,
                       max_running_time=None,
                       max_memory_per_node=None,
                       max_memmory_per_cpu=None,
                       modules_list=None,
                       environment_variables_dict=None):

        common_options = self.parse_common_options(tree_file=tree_file, output_format=output_format, show_xml=show_xml,
                                                   show_tree=show_tree, show_ancestral_sequences=show_ancestral_sequences,
                                                   show_evolutionary_events=show_evolutionary_events, showall=showall,
                                                   compute_posterior_support=compute_posterior_support, njtree=njtree,
                                                   skip_insertions=skip_insertions, codon_alignment=codon_alignment,
                                                   translated_alignment=translated_alignment)

        FileRoutines.safe_mkdir(output_directory)
        options_list = []
        for filename in list_of_files:
            basename = FileRoutines.split_filename(filename)[1]
            op = common_options
            op += " -d=%s" % filename
            op += " -o=%s/%s.fasta" % (output_directory,
                                       ("%s_%s" % (basename, output_suffix)) if output_suffix else basename)
            options_list.append(op)
        if handling_mode == "local":
            self.parallel_execute(options_list)
        elif handling_mode == "slurm":

            cmd_list = ["%s%s %s" % ((self.path + "/") if self.path else "", self.cmd, options) for options in options_list]
            self.slurm_run_multiple_jobs_in_wrap_mode(cmd_list,
                                                      cmd_log_file,
                                                      max_jobs=max_jobs,
                                                      job_name=job_name,
                                                      log_prefix=log_prefix,
                                                      error_log_prefix=error_log_prefix,
                                                      cpus_per_node=None,
                                                      max_running_jobs=None,
                                                      max_running_time=max_running_time,
                                                      cpus_per_task=cpus_per_task,
                                                      max_memory_per_node=max_memory_per_node,
                                                      max_memmory_per_cpu=max_memmory_per_cpu,
                                                      modules_list=modules_list,
                                                      environment_variables_dict=environment_variables_dict)

    def parallel_codon_alignment(self, list_of_files, output_directory, output_suffix=None, tree_file=None, output_format=None, show_xml=None,
                                 show_tree=None, show_ancestral_sequences=None, show_evolutionary_events=None,
                                 showall=None, compute_posterior_support=None, njtree=None,
                                 cmd_log_file=None,
                                 cpus_per_task=1,
                                 handling_mode="local",
                                 job_name=None,
                                 log_prefix=None,
                                 error_log_prefix=None,
                                 max_jobs=None,
                                 max_running_time=None,
                                 max_memory_per_node=None,
                                 max_memmory_per_cpu=None,
                                 modules_list=None,
                                 environment_variables_dict=None):

        self.parallel_align(list_of_files, output_directory, output_suffix=output_suffix, tree_file=tree_file,
                            output_format=output_format, show_xml=show_xml, show_tree=show_tree,
                            show_ancestral_sequences=show_ancestral_sequences,
                            show_evolutionary_events=show_evolutionary_events, showall=showall,
                            compute_posterior_support=compute_posterior_support, njtree=njtree, skip_insertions=True,
                            codon_alignment=True,
                            cmd_log_file=cmd_log_file,
                            cpus_per_task=cpus_per_task,
                            handling_mode=handling_mode,
                            job_name=job_name,
                            log_prefix=log_prefix,
                            error_log_prefix=error_log_prefix,
                            max_jobs=max_jobs,
                            max_running_time=max_running_time,
                            max_memory_per_node=max_memory_per_node,
                            max_memmory_per_cpu=max_memmory_per_cpu,
                            modules_list=modules_list,
                            environment_variables_dict=environment_variables_dict)
