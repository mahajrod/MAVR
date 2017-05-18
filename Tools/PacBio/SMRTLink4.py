#!/usr/bin/env python

from Tools.Abstract import Tool


class CCS(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "ccs", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(max_subread_len=None, min_subread_length=None, min_subread_num_per_circle=None):

        options = ""
        options += " --maxLength %i" % max_subread_len if max_subread_len else ""
        options += " --minLength %i" % min_subread_length if min_subread_length else ""
        options += " --minPasses %i" % min_subread_num_per_circle if min_subread_num_per_circle else ""

        return options

    def get_consensus(self, input_file, output_file, max_subread_len=None, min_subread_length=None,
                      min_subread_num_per_circle=None):

        options = self.parse_common_options(max_subread_len=max_subread_len,
                                            min_subread_length=min_subread_length,
                                            min_subread_num_per_circle=min_subread_num_per_circle)
        options += " %s" % input_file
        options += " %s" % output_file

        self.execute(options=options)
