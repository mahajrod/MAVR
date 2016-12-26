#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
from collections import Iterable


from Routines import SequenceRoutines

print_mutex = mp.Lock()


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use mutex to safe write to stdout from multiple threads
    print_mutex.acquire()
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    print_mutex.release()

    os.system(exe_string)


class Tool(SequenceRoutines):

    def __init__(self, cmd, path="", max_threads=4, jar_path=None, jar=None,
                 max_memory="500m", timelog="tool_time.log"):
        SequenceRoutines.__init__(self)
        self.path = self.check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        #print(jar_path)
        self.jar_path = self.check_path(jar_path) if jar_path else None
        self.jar = jar
        self.max_memory = max_memory
        self.timelog = timelog

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
                         write_output_to_file=None, external_process_pool=None, async_run=False):
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

        results = process_pool.map_async(execute, exe_string_list) if async_run else process_pool.map(execute, exe_string_list)
        #process_pool.close()
        if write_output_to_file:
            with open(write_output_to_file, "w") as out_fd:
                out_fd.write(results)
        return results if capture_output else None

    
class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):

        Tool.__init__(self, "java", path=java_path, max_threads=max_threads,
                      jar_path=jar_path, jar=jar, max_memory=max_memory, timelog=timelog)

    def execute(self, options="", cmd=None, capture_output=False):
        command = cmd if cmd is not None else ""

        java_string = "java"
        java_string += " -Xmx%s" % str(self.max_memory) if self.max_memory else ""
        #print (self.jar_path)
        java_string += " -jar %s%s" % (self.jar_path, self.jar)
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

