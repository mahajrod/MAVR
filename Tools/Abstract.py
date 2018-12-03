#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
from collections import Iterable

import numpy as np

from Routines.Sequence import SequenceRoutines
from Routines.Alignment import AlignmentRoutines

print_mutex = mp.Lock()


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use mutex to safe write to stdout from multiple threads
    print_mutex.acquire()
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    print_mutex.release()

    os.system(exe_string)


class Tool(SequenceRoutines, AlignmentRoutines):

    def __init__(self, cmd, path="", max_threads=4, jar_path="", jar=None,
                 max_memory="500m", max_per_thread_memory="500m", timelog=None, tmp_dir=None):
        SequenceRoutines.__init__(self)
        self.path = self.check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        #print(jar_path)
        self.jar_path = self.check_path(jar_path) if jar_path else ""
        self.jar = jar
        self.max_memory = max_memory
        self.timelog = timelog
        self.max_per_thread_memory = max_per_thread_memory
        self.tmp_dir = tmp_dir

    def execute(self, options="", cmd=None, directory=None, capture_output=False, generate_cmd_string_only=False, intepreter=None):
        command = cmd if cmd is not None else self.cmd

        exe_string = ""
        exe_string += " cd %s && " % self.check_dir_path(directory) if directory else ""
        exe_string += ("%s " % intepreter if intepreter else "") + (self.check_path(self.path) if self.path else "") + command + " " + options

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        exe_string = "time -f 'Time\\t%%E real,\\t%%U user,\\t%%S sys' -a -o %s %s" % (self.timelog, exe_string) if self.timelog else exe_string

        if generate_cmd_string_only:
            return exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None, dir_list=None,
                         write_output_to_file=None, external_process_pool=None, async_run=False,
                         job_chunksize=1, output_file_list=None, duplicate_to_stdout=False, intepreter=None):

        command = cmd if cmd is not None else self.cmd
        number_of_commands = len(options_list)

        if dir_list:
            if isinstance(dir_list, str):
                directory_list = [dir_list for i in range(0, number_of_commands)]
            elif len(dir_list) == number_of_commands:
                directory_list = [directory for directory in dir_list]
            else:
                raise ValueError("Error during option parsing for parallel execution in different folders. "
                                 "Length of directory list is not equal to length of option list")
        else:
            directory_list = [None for i in range(0, number_of_commands)]

        if output_file_list:
            if len(dir_list) == number_of_commands:
                out_list = [out_file for out_file in output_file_list]
            else:
                raise ValueError("Error during option parsing for parallel execution with storing stdout and stderr. "
                                 "Length of output file list is not equal to length of option list")
        else:
            out_list = [None for i in range(0, number_of_commands)]

        exe_string_list = []

        for options, directory, output in zip(options_list, directory_list, out_list):
            com = ""
            com += " cd %s && " % self.check_dir_path(directory) if directory else ""
            com += ("%s " % intepreter if intepreter else "") + " %s%s" % (self.check_path(self.path) if self.path else "", command)
            com += " %s" % options

            if output:
                com += " | tee %s/%s" % (directory, output) if duplicate_to_stdout else " > %s/%s 2>&1" % (directory, output)

            exe_string_list.append(com)
        """
        if dir_list:
            if isinstance(dir_list, str):
                if output_file_list:
                    for options, output in zip(options_list, output_file_list):
                        com = ""
                        com += " cd %s && " % dir_list
                        com += " %s" % self.check_path(self.path) if self.path else ""
                        com += " %s" % command
                        com += " %s" % options
                        com += " | tee %s/%s" % (dir_list, output) if duplicate_to_stdout else " > %s/%s 2>&1" % (dir_list, output)

                        exe_string_list.append(com)

                        #exe_string_list = [("cd %s && " % dir_list) + (self.check_path(self.path) if self.path else "")
                        #               + command + " " + options +
                        #               (" | tee %s/%s" % (dir_list, output) if duplicate_to_stdout else " > %s/%s 2>&1" % (dir_list, output)) for options, output in zip(options_list, output_file_list)]
                else:
                    for options in options_list:
                        com = ""
                        com += " cd %s && " % dir_list
                        com += " %s" % self.check_path(self.path) if self.path else ""
                        com += " %s" % command
                        com += " %s" % options

                        exe_string_list.append(com)

                        #exe_string_list = [("cd %s && " % dir_list) + (self.check_path(self.path) if self.path else "")
                        #               + command + " " + options for options in options_list]

            elif isinstance(dir_list, Iterable) and (len(options_list) == len(dir_list)):
                if output_file_list:
                    for options, directory, output in zip(options_list, dir_list, output_file_list):
                        com = ""
                        com += " cd %s && " % directory
                        com += " %s" % self.check_path(self.path) if self.path else ""
                        com += " %s" % command
                        com += " %s" % options
                        com += " | tee %s/%s" % (directory, output) if duplicate_to_stdout else " > %s/%s 2>&1" % (directory, output)

                        exe_string_list.append(com)

                        #exe_string_list = [("cd %s && " % directory) + (self.check_path(self.path) if self.path else "")
                        #               + command + " " + options +
                        #               (" | tee %s/%s" % (directory, output) if duplicate_to_stdout else " > %s/%s 2>&1" % (directory, output)) for options, directory, output in zip(options_list, dir_list, output_file_list)]

                else:
                    for options, directory in zip(options_list, dir_list):
                        com = ""
                        com += " cd %s && " % directory
                        com += " %s" % self.check_path(self.path) if self.path else ""
                        com += " %s" % command
                        com += " %s" % options

                        exe_string_list.append(com)

                        #exe_string_list = [("cd %s && " % directory) + (self.check_path(self.path) if self.path else "")
                        #               + command + " " + options for options, directory in zip(options_list, dir_list)]

            else:
                raise ValueError("Error during option parsing for parallel execution in different folders. "
                                 "Length of directory list is not equal to length of option list")
        elif output_file_list:
            for options, output in zip(options_list, output_file_list):
                com = ""
                com += " %s" % self.check_path(self.path) if self.path else ""
                com += " %s" % command
                com += " %s" % options
                com += " | tee %s" % output if duplicate_to_stdout else " > %s 2>&1" % output

                exe_string_list.append(com)

                #exe_string_list = [(self.check_path(self.path) if self.path else "") + command + " " + options +
                #               (" | tee %s" % output if duplicate_to_stdout else " > %s 2>&1" % output) for options, output in zip(options_list, output_file_list)]

        else:
            for options in options_list:
                com = ""
                com += " %s" % self.check_path(self.path) if self.path else ""
                com += " %s" % command
                com += " %s" % options

                exe_string_list.append(com)

                #exe_string_list = [(self.check_path(self.path) if self.path else "") + command + " " + options
                #               for options in options_list]

        """

        with open("exe_list.t", "a") as exe_fd:
            for entry in exe_string_list:
                exe_fd.write("%s\n" % entry)
        process_pool = external_process_pool if external_process_pool else mp.Pool(threads if threads else self.threads)

        results = process_pool.map_async(execute, exe_string_list, chunksize=job_chunksize) if async_run else process_pool.map(execute, exe_string_list, chunksize=job_chunksize)
        #process_pool.close()
        if write_output_to_file:
            with open(write_output_to_file, "w") as out_fd:
                out_fd.write(results)
        return results if capture_output else None

    def generate_slurm_job_script(self,
                                  job_name,
                                  log_prefix,
                                  task_commands,
                                  error_log_prefix,
                                  node_number=None,
                                  partition=None,
                                  qos=None,
                                  cpus_per_node=None,
                                  node_type=None,
                                  max_memory_per_node=None,
                                  stdout_file=None,
                                  email=None,
                                  mail_type=None,
                                  job_script=None,
                                  task_index_list=None,
                                  start_task_index=None,
                                  end_task_index=None,
                                  max_running_jobs=None,
                                  max_running_time=None,
                                  cpus_per_task=None,
                                  max_memmory_per_cpu=None,
                                  modules_list=None,
                                  environment_variables_dict=None):

        #if (not task_index_list) and ((not start_task_index) or (not end_task_index)):
        #    raise ValueError("Neither task index list nor start or end task index were set")

        script = "#!/usr/bin/env bash\n"
        if task_index_list or start_task_index or end_task_index:
            script += "#SBATCH --array=%s" % (",".join(map(str, task_index_list)) if task_index_list else "%s-%s" % (str(start_task_index),
                                                                                                                     str(end_task_index)))
            script += "%s\n" % ("%%%i" % max_running_jobs if max_running_jobs else "")

        script += "#SBATCH --time=%s         # Run time in hh:mm:ss\n" % max_running_time if max_running_time else ""
        script += "#SBATCH --mem-per-cpu=%i       # Minimum memory required per CPU (in megabytes)\n" % max_memmory_per_cpu if max_memmory_per_cpu else ""
        script += "#SBATCH --cpus-per-task=%i\n" % cpus_per_task if cpus_per_task else ""
        script += "#SBATCH --job-name=%s\n" % job_name if job_name else ""
        script += "#SBATCH --error=%s.%%A_%%a.err\n" % error_log_prefix
        script += "#SBATCH --output=%s.%%A_%%a.log\n" % log_prefix

        script += "#SBATCH --nodes=%i\n" % node_number if node_number else ""
        script += "#SBATCH --partition=%s\n" % partition if partition else ""
        script += "#SBATCH --ntasks-per-nod=%i\n" % cpus_per_node if cpus_per_node else ""

        script += "#SBATCH --node_type=%s\n" % node_type if node_type else ""

        script += "#SBATCH --mem=%s\n" % str(max_memory_per_node) if max_memory_per_node else ""
        script += "#SBATCH --mail-user=%s\n" % email if email else ""
        script += "#SBATCH --mail-type=%s\n" % mail_type if mail_type else ""
        script += "#SBATCH --output=%s\n" % stdout_file if stdout_file else ""

        script += "\n"

        if environment_variables_dict:
            script += "#--------------Environment variables--------------\n\n"
            for variable in environment_variables_dict:
                script += "%s=%s\n" % (variable, environment_variables_dict[variable])
            script += "#-------------------------------------------\n\n"

        if modules_list:
            script += "#--------------Modules--------------\n\n"
            modules_to_load = [modules_list] if isinstance(modules_list, str) else list(modules_list)
            for module in modules_to_load:
                script += "module load %s\n" % module

            script += "#-------------------------------------------\n\n"

        script += "#--------------Commands--------------\n\n"

        script += "%s\n" % task_commands

        script += "#-------------------------------------------\n\n"

        if job_script:
            with open(job_script, "w") as job_fd:
                job_fd.write(script)

        return script

    def slurm_run_job_from_script(self,
                                  job_array_script,
                                  after_job_id_list=[],
                                  afterany_job_id_list=[],
                                  afternotok_job_id_list=[],
                                  afterok_job_id_list=[],
                                  singleton=False):
        """
        after:jobid[:jobid...]	    job can begin after the specified jobs have started
        afterany:jobid[:jobid...]	job can begin after the specified jobs have terminated
        afternotok:jobid[:jobid...]	job can begin after the specified jobs have failed
        afterok:jobid[:jobid...]	job can begin after the specified jobs have run to completion with an exit code of zero (see the user guide for caveats).
        singleton	                jobs can begin execution after all previously launched jobs with the same name and user have ended. This is useful to collate results of a swarm or to send a notification at the end of a swarm.

        """

        dependencies_option_list = []

        if after_job_id_list:
            dependencies_option_list.append("after:%s" % (after_job_id_list if isinstance(after_job_id_list, str) else ":".join(after_job_id_list)))
        if afterany_job_id_list:
            dependencies_option_list.append("afterany:%s" % (afterany_job_id_list if isinstance(afterany_job_id_list, str) else ":".join(afterany_job_id_list)))
        if afternotok_job_id_list:
            dependencies_option_list.append("afternotok:%s" % (afternotok_job_id_list if isinstance(afternotok_job_id_list, str) else ":".join(afternotok_job_id_list)))
        if afterok_job_id_list:
            #print afterok_job_id_list
            dependencies_option_list.append("afterok:%s" % (afterok_job_id_list if isinstance(afterok_job_id_list, str) else ":".join(afterok_job_id_list)))
        if singleton:
            dependencies_option_list.append("singleton")

        options = ""

        if dependencies_option_list:
            dependencies_option = ",".join(dependencies_option_list)
            options += " --dependency=%s" % dependencies_option

        options += " %s" % job_array_script

        command = "sbatch %s " % options

        # Popen.stdout returns file object
        print(command)
        job_id = Popen([command], shell=True, stdout=PIPE).stdout.readline().strip().split()[-1]
        #print "Job id", job_id
        return job_id

    def slurm_run_job(self,
                      job_name,
                      log_prefix,
                      task_commands,
                      error_log_prefix,
                      job_script,
                      node_number=None,
                      partition=None,
                      qos=None,
                      cpus_per_node=None,
                      node_type=None,
                      max_memory_per_node=None,
                      stdout_file=None,
                      email=None,
                      mail_type=None,
                      task_index_list=None,
                      start_task_index=None,
                      end_task_index=None,
                      max_running_jobs=None,
                      max_running_time=None,
                      cpus_per_task=None,
                      max_memmory_per_cpu=None,
                      modules_list=None,
                      environment_variables_dict=None,
                      after_job_id_list=[],
                      afterany_job_id_list=[],
                      afternotok_job_id_list=[],
                      afterok_job_id_list=[],
                      singleton=False):

        self.generate_slurm_job_script(job_name,
                                       log_prefix,
                                       task_commands,
                                       error_log_prefix,
                                       node_number=node_number,
                                       partition=partition,
                                       qos=qos,
                                       cpus_per_node=cpus_per_node,
                                       node_type=node_type,
                                       max_memory_per_node=max_memory_per_node,
                                       stdout_file=stdout_file,
                                       email=email,
                                       mail_type=mail_type,
                                       job_script=job_script,
                                       task_index_list=task_index_list,
                                       start_task_index=start_task_index,
                                       end_task_index=end_task_index,
                                       max_running_jobs=max_running_jobs,
                                       max_running_time=max_running_time,
                                       cpus_per_task=cpus_per_task,
                                       max_memmory_per_cpu=max_memmory_per_cpu,
                                       modules_list=modules_list,
                                       environment_variables_dict=environment_variables_dict)

        return self.slurm_run_job_from_script(job_script,
                                              after_job_id_list=after_job_id_list,
                                              afterany_job_id_list=afterany_job_id_list,
                                              afternotok_job_id_list=afternotok_job_id_list,
                                              afterok_job_id_list=afterok_job_id_list,
                                              singleton=singleton)

    def parse_sbatch_options(self,
                             job_name=None,
                             log_prefix=None,
                             error_log_prefix=None,
                             node_number=None,
                             partition=None,
                             cpus_per_node=None,
                             node_type=None,
                             max_memory_per_node=None,
                             stdout_file=None,
                             email=None,
                             mail_type=None,
                             task_index_list=None,
                             start_task_index=None,
                             end_task_index=None,
                             max_running_jobs=None,
                             max_running_time=None,
                             cpus_per_task=None,
                             max_memmory_per_cpu=None):

        options = ""
        if task_index_list or start_task_index or end_task_index:
            options += " --array=%s" % (",".join(map(str, task_index_list)) if task_index_list else "%s-%s" % (str(start_task_index),
                                                                                                                     str(end_task_index)))
            options += "%s" % ("%%%i" % max_running_jobs if max_running_jobs else "")

        options += " --time=%s" % max_running_time if max_running_time else ""
        options += " --mem-per-cpu=%i" % max_memmory_per_cpu if max_memmory_per_cpu else ""
        options += " --cpus-per-task=%i" % cpus_per_task if cpus_per_task else ""
        options += " --job-name=%s" % job_name if job_name else ""
        options += " --error=%s.%%A_%%a.err" % error_log_prefix if error_log_prefix else ""
        options += " --output=%s.%%A_%%a.log" % log_prefix if log_prefix else ""

        options += " --nodes=%i" % node_number if node_number else ""
        options += " --partition=%s" % partition if partition else ""
        options += " --ntasks-per-nod=%i" % cpus_per_node if cpus_per_node else ""

        options += " --node_type=%s" % node_type if node_type else ""

        options += " --mem=%s" % str(max_memory_per_node) if max_memory_per_node else ""
        options += " --mail-user=%s" % email if email else ""
        options += " --mail-type=%s" % mail_type if mail_type else ""
        options += " --output=%s" % stdout_file if stdout_file else ""

        return options

    def slurm_run_multiple_jobs_in_wrap_mode(self,
                                             cmd_list,
                                             cmd_log,
                                             max_jobs=None,
                                             job_name=None,
                                             log_prefix=None,
                                             error_log_prefix=None,
                                             node_number=None,
                                             partition=None,
                                             cpus_per_node=None,
                                             node_type=None,
                                             max_memory_per_node=None,
                                             stdout_file=None,
                                             email=None,
                                             mail_type=None,
                                             task_index_list=None,
                                             start_task_index=None,
                                             end_task_index=None,
                                             max_running_jobs=None,
                                             max_running_time=None,
                                             cpus_per_task=None,
                                             max_memmory_per_cpu=None,
                                             max_memory_per_cpu_per_task_list=None,
                                             modules_list=None,
                                             environment_variables_dict=None):
        if max_memory_per_cpu_per_task_list is not None:
            if len(max_memory_per_cpu_per_task_list) != cmd_list:
                raise ValueError("Length of cmd_list(%i) is noq equal to max_memory_per_cpu_per_task_list(%i)" % (len(cmd_list), len(max_memory_per_cpu_per_task_list)))

        with open(cmd_log, "w") as cmd_fd:
            cmd_fd.write("#job_id\tcmd_ids\tsbatch_cmd\n")
            environment_variables_cmd = ""
            modules_cmd = ""

            if environment_variables_dict:
                for variable in environment_variables_dict:
                    environment_variables_cmd += " %s=%s; " % (variable, environment_variables_dict[variable])

            if modules_list:
                modules_to_load = [modules_list] if isinstance(modules_list, str) else list(modules_list)
                for module in modules_to_load:
                    modules_cmd += " module load %s; " % module

            sbatch_options = self.parse_sbatch_options(job_name=job_name,
                                                       log_prefix=log_prefix,
                                                       error_log_prefix=error_log_prefix,
                                                       node_number=node_number,
                                                       partition=partition,
                                                       cpus_per_node=cpus_per_node,
                                                       node_type=node_type,
                                                       max_memory_per_node=max_memory_per_node,
                                                       stdout_file=stdout_file,
                                                       email=email,
                                                       mail_type=mail_type,
                                                       task_index_list=task_index_list,
                                                       start_task_index=start_task_index,
                                                       end_task_index=end_task_index,
                                                       max_running_jobs=max_running_jobs,
                                                       max_running_time=max_running_time,
                                                       cpus_per_task=cpus_per_task,
                                                       max_memmory_per_cpu=max_memmory_per_cpu)

            job_submit_index_array = np.linspace(0, len(cmd_list), max_jobs, dtype=int) if max_jobs else [i for i in range(0, len(cmd_list) + 1)]

            for job_index in range(0, len(job_submit_index_array) - 1):
                if max_memory_per_cpu_per_task_list:
                    max_memory_per_cpu_for_job = max(max_memory_per_cpu_per_task_list[job_submit_index_array[job_index]:job_submit_index_array[job_index+1]])
                    sbatch_options = self.parse_sbatch_options(job_name=job_name,
                                                               log_prefix=log_prefix,
                                                               error_log_prefix=error_log_prefix,
                                                               node_number=node_number,
                                                               partition=partition,
                                                               cpus_per_node=cpus_per_node,
                                                               node_type=node_type,
                                                               max_memory_per_node=max_memory_per_node,
                                                               stdout_file=stdout_file,
                                                               email=email,
                                                               mail_type=mail_type,
                                                               task_index_list=task_index_list,
                                                               start_task_index=start_task_index,
                                                               end_task_index=end_task_index,
                                                               max_running_jobs=max_running_jobs,
                                                               max_running_time=max_running_time,
                                                               cpus_per_task=cpus_per_task,
                                                               max_memmory_per_cpu=max_memory_per_cpu_for_job)

                options = "%s" % sbatch_options
                options += " --wrap \"%s %s %s\"" % (environment_variables_cmd,
                                                     modules_cmd,
                                                     ";".join(cmd_list[job_submit_index_array[job_index]:job_submit_index_array[job_index+1]]))
                command = "sbatch %s" % options

                print command

                job_id = Popen([command], shell=True, stdout=PIPE).stdout.readline().strip().split()[-1]
                cmd_fd.write("%s\t%i-%i\t%s\n" % (job_id,
                                                  job_submit_index_array[job_index],
                                                  job_submit_index_array[job_index+1] - 1, command))


class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):

        Tool.__init__(self, "java", path=java_path, max_threads=max_threads,
                      jar_path=jar_path, jar=jar, max_memory=max_memory, timelog=timelog, tmp_dir=None)

    def execute(self, options="", cmd=None, capture_output=False, runtype="jar", generate_cmd_string_only=False):
        command = cmd if cmd is not None else ""

        java_string = "java"
        java_string += " -Xmx%s" % str(self.max_memory) if self.max_memory else ""
        java_string += " -Djava.io.tmpdir=%s" % self.tmp_dir if self.tmp_dir else ""
        #print (self.jar_path)
        java_string += " -%s %s%s" % (runtype, self.check_dir_path(self.jar_path) if self.jar_path else "", self.jar)
        java_string += " %s" % command
        java_string += " %s" % options

        exe_string = (self.check_path(self.path) if self.path else "") + java_string

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        exe_string = "time -f 'Time\\t%%E real,\\t%%U user,\\t%%S sys' -a -o %s %s" % (self.timelog, exe_string) if self.timelog else exe_string

        if generate_cmd_string_only:
            return exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None, dir_list=None,
                         write_output_to_file=None, external_process_pool=None, async_run=False, runtype="jar",
                         job_chunksize=1):
        command = cmd if cmd is not None else self.cmd
        if dir_list:
            if isinstance(dir_list, str):
                exe_string_list = [("cd %s && " % dir_list) + (self.check_path(self.path) if self.path else "")
                                   + command + " -%s %s/%s" % (runtype, self.check_dir_path(self.jar_path) if self.jar_path else "", self.jar)
                                   + " " + options for options in options_list]

            elif isinstance(dir_list, Iterable) and (len(options_list) == len(dir_list)):
                exe_string_list = [("cd %s && " % directory) + (self.check_path(self.path) if self.path else "")
                                   + command + " -%s %s/%s" % (runtype, self.check_dir_path(self.jar_path) if self.jar_path else "", self.jar)
                                   + " " + options for options, directory in zip(options_list, dir_list)]
            else:
                raise ValueError("Error during option parsing for parallel execution in different folders. "
                                 "Length of directory list is not equal to length of option list")

        else:
            exe_string_list = [(self.check_path(self.path) if self.path else "")
                               + command + " -%s %s/%s" % (runtype, self.check_dir_path(self.jar_path) if self.jar_path else "", self.jar)
                               + " " + options for options in options_list]

        with open("exe_list.t", "a") as exe_fd:
            for entry in exe_string_list:
                exe_fd.write("%s\n" % entry)
        process_pool = external_process_pool if external_process_pool else mp.Pool(threads if threads else self.threads)

        results = process_pool.map_async(execute, exe_string_list, chunksize=job_chunksize) if async_run else process_pool.map(execute, exe_string_list, chunksize=job_chunksize)
        #process_pool.close()
        if write_output_to_file:
            with open(write_output_to_file, "w") as out_fd:
                out_fd.write(results)
        return results if capture_output else None