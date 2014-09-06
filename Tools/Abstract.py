#!/usr/bin/env python

import os
from Routines.Functions import check_path


class Tool():

    def __init__(self, cmd, path="", max_threads=4, jar_path=None, jar=None, max_memory="500m"):
        self.path = check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        self.jar_path = jar_path
        self.jar = jar
        self.max_memory = max_memory

    def execute(self, options, cmd=None):
        command = cmd if cmd is not None else self.cmd
        if self.jar_path:
            exe_string = self.path + "java -Xmx%s -jar %s%s %s" % (self.max_memory, self.jar_path, self.jar, command) \
                         + " " + options
        else:
            exe_string = self.path + command + " " + options
        print("Executing:\n\t%s" % exe_string)
        os.system(exe_string)