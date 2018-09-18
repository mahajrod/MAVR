#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
from collections import Iterable


from Routines.Sequence import SequenceRoutines

print_mutex = mp.Lock()


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use mutex to safe write to stdout from multiple threads
    print_mutex.acquire()
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    print_mutex.release()

    os.system(exe_string)


class Tool(SequenceRoutines):

    def __init__(self, cmd, path="", max_threads=4, jar_path="", jar=None,
                 max_memory="500m", max_per_thread_memory="500m", timelog=None):
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

    def execute(self, options="", cmd=None, capture_output=False):
        command = cmd if cmd is not None else self.cmd

        exe_string = (self.check_path(self.path) if self.path else "") + command + " " + options

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        exe_string = "time -f 'Time\\t%%E real,\\t%%U user,\\t%%S sys' -a -o %s %s" % (self.timelog, exe_string) if self.timelog else exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None, dir_list=None,
                         write_output_to_file=None, external_process_pool=None, async_run=False,
                         job_chunksize=1):
        command = cmd if cmd is not None else self.cmd
        if dir_list:
            if isinstance(dir_list, str):
                exe_string_list = [("cd %s && " % dir_list) + (self.check_path(self.path) if self.path else "")
                                   + command + " " + options for options in options_list]
            elif isinstance(dir_list, Iterable) and (len(options_list) == len(dir_list)):
                exe_string_list = [("cd %s && " % directory) + (self.check_path(self.path) if self.path else "")
                                   + command + " " + options for options, directory in zip(options_list, dir_list)]
            else:
                raise ValueError("Error during option parsing for parallel execution in different folders. "
                                 "Length of directory list is not equal to length of option list")

        else:
            exe_string_list = [(self.check_path(self.path) if self.path else "") + command + " " + options
                               for options in options_list]

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

    def generate_slurm_job_array_script(self,
                                        job_name,
                                        log_prefix,
                                        task_commands,
                                        error_log_prefix,
                                        job_array_script_file=None,
                                        task_index_list=None,
                                        start_task_index=None,
                                        end_task_index=None,
                                        max_running_jobs=None,
                                        max_running_time=None,
                                        max_memmory_per_cpu=None):

        if (not task_index_list) and ((not start_task_index) or (not end_task_index)):
            raise ValueError("Neither task index list nor start or end task index were set")

        script = "#!/usr/bin/env bash\n"

        script += "#SBATCH --array=%s" % (",".join(map(str, task_index_list)) if task_index_list else "%s-%s" % (str(start_task_index),
                                                                                                                 str(end_task_index)))
        script += "%s\n" % ("%%%i" % max_running_jobs if max_running_jobs else "")
        script += "#SBATCH --time=%s         # Run time in hh:mm:ss\n" % max_running_time if max_running_time else ""
        script += "#SBATCH --mem-per-cpu=%i       # Minimum memory required per CPU (in megabytes)\n" % max_memmory_per_cpu if max_memmory_per_cpu else ""
        script += "#SBATCH --job-name=%s\n" % job_name if job_name else ""
        script += "#SBATCH --error=%s.%%A_%%a.err\n" % error_log_prefix
        script += "#SBATCH --output=%s.%%A_%%a.log\n" % log_prefix

        script += "%s\n" % task_commands

        if job_array_script_file:
            with open(job_array_script_file, "w") as job_fd:
                job_fd.write(script)

        return script

    def slurm_run_job_array(self, job_array_script):

        # Popen.stdout returns file object
        job_id = Popen(["sbatch %s" % job_array_script], shell=True, stdout=PIPE).stdout.readline().strip().split()[-1]

        return job_id


class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):

        Tool.__init__(self, "java", path=java_path, max_threads=max_threads,
                      jar_path=jar_path, jar=jar, max_memory=max_memory, timelog=timelog)

    def execute(self, options="", cmd=None, capture_output=False, runtype="jar"):
        command = cmd if cmd is not None else ""

        java_string = "java"
        java_string += " -Xmx%s" % str(self.max_memory) if self.max_memory else ""
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