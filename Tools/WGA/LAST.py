#!/usr/bin/env python

import os
from Tools.Abstract import Tool
from Routines import DrawingRoutines


class LAST(Tool):
    """
    Class for samtools 1.0+
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "last", path=path, max_threads=max_threads)

    def parse_lastdb_options(self, output, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                             verbose=True, keep_preliminary_masking=True,
                             mask_simple_repeats=True):
        options = " -P %i" % self.threads
        options += " -c" if softmasking else ""
        options += " -u %s" % seeding_scheme if seeding_scheme else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)

        options += " -v" if verbose else verbose

        options += " %s" % output
        options += " %s" % (input_fasta_list if isinstance(input_fasta_list, str) else " ".join(input_fasta_list))

        return options

    def create_last_db(self, output, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                       verbose=True, keep_preliminary_masking=True, mask_simple_repeats=True):

        options = self.parse_lastdb_options(output, input_fasta_list=input_fasta_list,
                                            softmasking=softmasking,
                                            seeding_scheme=seeding_scheme, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats)

        self.execute(options=options, cmd="lastdb")

