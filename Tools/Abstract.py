#!/usr/bin/env python

import os
from Routines.Functions import check_path


class Tool():

    def __init__(self, cmd, path=""):
        self.path = check_path(path)
        self.cmd = cmd

    def execute(self, options, command=None):
        if command:
            cmd = command
        else:
            cmd = self.cmd
        exe_string = self.path + cmd + " " + options
        print("Executing:\n\t%s" % exe_string)
        os.system(exe_string)