#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Routines.File import save_mkdir, split_filename, check_path
from Tools.Abstract import Tool


class BioConvert(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "bfconvert", path=path, max_threads=max_threads)

    def convert(self, input_file, output_file):

        options = " %i" % input_file
        options += " %s" % output_file

        self.execute(options)

    def parallel_convert(self, list_of_files, output_directory, output_format=".tiff"):

        save_mkdir(output_directory)
        options_list = []

        for filename in list_of_files:
            option = " %s" % filename
            option += " %s%s%s" % (check_path(output_directory), split_filename(filename)[1], output_format)
            options_list.append(option)

        self.parallel_execute(options_list)

if __name__ == "__main__":
    pass
