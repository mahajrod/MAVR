#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class FastQC(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "fastqc", path=path, max_threads=max_threads)

    def check(self, input_files, out_dir="./", nogroup=True, noextract=False, adapters_file=None):
        options = " -o %s" % out_dir
        options += " --nogroup" if nogroup else ""
        options += " --noextract" if noextract else ""
        options += " -t %i" % self.threads
        options += " -a %s" % adapters_file if adapters_file else ""
        options += " %s" % input_files if isinstance(input_files, str) else " ".join(input_files)
        self.execute(options, cmd="fastqc")


if __name__ == "__main__":
    pass