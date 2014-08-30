#!/usr/bin/env python

import os
from Routines.Functions import check_path


class Tool():

    def __init__(self, cmd, path="", max_threads=4):
        self.path = check_path(path)
        self.cmd = cmd
        self.threads = max_threads

    def execute(self, options, cmd=None):
        command = cmd if cmd else self.cmd
        exe_string = self.path + command + " " + options
        print("Executing:\n\t%s" % exe_string)
        os.system(exe_string)