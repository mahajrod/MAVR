#!/usr/bin/env python
import os

from Routines import FileRoutines

from Tools.Abstract import Tool


class SWAMP(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "SWAMP.py", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(input_folder=None, branch_names_file=None, window_size=None,
                             substitution_threshold=None, min_seq_len=None, interscan_masking=False):

        options = " -i %s" % input_folder if input_folder else ""
        options += " -b %s" % branch_names_file if branch_names_file else ""
        options += " -t %i" % substitution_threshold if substitution_threshold else ""
        options += " -w %i" % window_size if window_size else ""
        options += " -m %i" % min_seq_len if min_seq_len else ""
        options += " -s" if interscan_masking else ""

    def mask_alignments(self, input_folder, output_folder, ):

        self.safe_mkdir(output_folder)
        for filename in os.listdir(input_folder):
            pass