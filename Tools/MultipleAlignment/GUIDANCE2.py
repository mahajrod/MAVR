#!/usr/bin/env python
import os
from Tools.Abstract import Tool

from Routines import FileRoutines


class GUIDANCE2(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "guidance.pl", path=path, max_threads=max_threads)

    def parse_common_options(self, seq_file=None, output_prefix=None, output_dir=None, msa_tool='prank',
                             seq_type=None, bootstrap_number=100, genetic_code=1, threads=None,
                             msa_tool_options=None, seq_cutoff=None, col_cutoff=None, mafft_bin=None,
                             prank_bin=None, muscle_bin=None, pagan_bin=None, ruby_bin=None, program=None):
        """
        1) Nuclear Standard
        15) Nuclear Blepharisma
        6) Nuclear Ciliate
        10) Nuclear Euplotid
        2) Mitochondria Vertebrate
        5) Mitochondria Invertebrate
        3) Mitochondria Yeast
        13) Mitochondria Ascidian
        9) Mitochondria Echinoderm
        14) Mitochondria Flatworm
        4) Mitochondria Protozoan
        """
        options = " --outDir %s" % output_dir if output_dir else ""
        options += " --program %s" % program if program else ""
        options += " --proc_num %i" % (threads if threads else self.threads)
        options += " --msaProgram %s" % msa_tool
        options += " --MSA_Param %s" % msa_tool_options if msa_tool_options else ""
        options += " --seqType %s" % seq_type
        options += " --bootstraps %i" % bootstrap_number if bootstrap_number else ""
        options += " --genCode %i" % genetic_code if genetic_code else ""
        options += " --seqCutoff %f" % seq_cutoff if seq_cutoff else ""
        options += " --colCutoff %f" % col_cutoff if col_cutoff else ""

        options += " --mafft %s" % mafft_bin if mafft_bin else ""
        options += " --prank %s" % prank_bin if prank_bin else ""
        options += " --muscle %s" % muscle_bin if muscle_bin else ""
        options += " --pagan %s" % pagan_bin if pagan_bin else ""
        options += " --ruby %s" % ruby_bin if ruby_bin else ""

        options += " --seqFile %s" % seq_file if seq_file else ""
        options += " --dataset %s" % output_prefix if output_prefix else ""

        return options

    def align(self, seq_file, output_prefix, output_dir, msa_tool='prank',
              seq_type=None, bootstrap_number=100, genetic_code=1, threads=None,
              msa_tool_options=None, seq_cutoff=None, col_cutoff=None, mafft_bin=None,
              prank_bin=None, muscle_bin=None, pagan_bin=None, ruby_bin=None, program=None):

        options = self.parse_common_options(seq_file=seq_file, output_prefix=output_prefix,
                                            output_dir=output_dir, msa_tool=msa_tool,
                                            seq_type=seq_type, bootstrap_number=bootstrap_number,
                                            genetic_code=genetic_code, threads=threads,
                                            msa_tool_options=msa_tool_options, seq_cutoff=seq_cutoff,
                                            col_cutoff=col_cutoff, mafft_bin=mafft_bin,
                                            prank_bin=prank_bin, muscle_bin=muscle_bin,
                                            pagan_bin=pagan_bin, ruby_bin=ruby_bin, program=program)
        self.execute(options)

    def parallel_align(self, list_of_files, output_dir, msa_tool='prank',
                       seq_type=None, bootstrap_number=100, genetic_code=1, threads=None,
                       msa_tool_options=None, seq_cutoff=None, col_cutoff=None, mafft_bin=None,
                       prank_bin=None, muscle_bin=None, pagan_bin=None, ruby_bin=None, program=None,
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

        common_options = self.parse_common_options(output_dir=output_dir, msa_tool=msa_tool,
                                                   seq_type=seq_type, bootstrap_number=bootstrap_number,
                                                   genetic_code=genetic_code, threads=threads,
                                                   msa_tool_options=msa_tool_options, seq_cutoff=seq_cutoff,
                                                   col_cutoff=col_cutoff, mafft_bin=mafft_bin,
                                                   prank_bin=prank_bin, muscle_bin=muscle_bin,
                                                   pagan_bin=pagan_bin, ruby_bin=ruby_bin, program=program)

        FileRoutines.safe_mkdir(output_dir)
        options_list = []
        for filename in list_of_files:
            basename = FileRoutines.split_filename(filename)[1]
            op = common_options
            op += " --seqFile %s" % filename
            op += " --dataset %s" % basename
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

    def parallel_codon_alignment(self, list_of_files, output_dir, msa_tool='prank',
                                 bootstrap_number=100, genetic_code=1, threads=None,
                                 msa_tool_options=None, seq_cutoff=None, col_cutoff=None, mafft_bin=None,
                                 prank_bin=None, muscle_bin=None, pagan_bin=None, ruby_bin=None, program=None,
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

        self.parallel_align(list_of_files, output_dir, msa_tool=msa_tool,
                            seq_type="codon", bootstrap_number=bootstrap_number,
                            genetic_code=genetic_code, threads=threads,
                            msa_tool_options=msa_tool_options, seq_cutoff=seq_cutoff, col_cutoff=col_cutoff,
                            mafft_bin=mafft_bin,
                            prank_bin=prank_bin, muscle_bin=muscle_bin,
                            pagan_bin=pagan_bin, ruby_bin=ruby_bin, program=program,
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
