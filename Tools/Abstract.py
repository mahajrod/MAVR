#!/usr/bin/env python

import os
from Routines.Functions import check_path


class Tool():

    def __init__(self, cmd, path=None, max_threads=4, jar_path=None, jar=None, max_memory="500m"):
        self.path = check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        self.jar_path = check_path(jar_path) if jar_path else None
        self.jar = jar
        self.max_memory = max_memory

    def execute(self, options="", cmd=None):

        command = cmd if cmd is not None else self.cmd
        exe_string = check_path(self.path) + command + " " + options
        print("Executing:\n\t%s" % exe_string)

        os.system(exe_string)


class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory="1g"):

        Tool.__init__(self, "java", path=check_path(java_path), max_threads=max_threads,
                      jar_path=check_path(jar_path),
                      jar=jar, max_memory=max_memory)

    def execute(self, options="", cmd=None):
        command = cmd if cmd is not None else ""
        exe_string = check_path(self.path) + "java -Xmx%s -jar %s%s %s %s" % (self.max_memory,
                                                                              check_path(self.jar_path),
                                                                              self.jar, command, options)
        print("Executing:\n\t%s" % exe_string)

        os.system(exe_string)

