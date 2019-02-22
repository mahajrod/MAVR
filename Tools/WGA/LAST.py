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

    def parse_lastdb_options(self, db_prefix, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                             verbose=True, keep_preliminary_masking=True,
                             mask_simple_repeats=True):
        options = " -P %i" % self.threads
        options += " -c" if softmasking else ""
        options += " -u %s" % seeding_scheme if seeding_scheme else ""
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)

        options += " -v" if verbose else verbose

        options += " %s" % db_prefix
        options += " %s" % (input_fasta_list if isinstance(input_fasta_list, str) else " ".join(input_fasta_list))

        return options

    def create_last_db(self, db_prefix, input_fasta_list, softmasking=True, seeding_scheme="YASS",
                       verbose=True, keep_preliminary_masking=True, mask_simple_repeats=True):

        options = self.parse_lastdb_options(db_prefix, input_fasta_list=input_fasta_list,
                                            softmasking=softmasking,
                                            seeding_scheme=seeding_scheme, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats)

        self.execute(options=options, cmd="lastdb")

    def parse_lastal_options(self, lastdb, query, output, verbose=True,
                             keep_preliminary_masking=True, mask_simple_repeats=True,
                             output_format="MAF", per_thread_memory="4G"):

        options = " -P %i" % self.threads

        options += " -v" if verbose else verbose
        options += " -R%i%i" % (1 if keep_preliminary_masking else 0,
                                1 if mask_simple_repeats else 0)
        options += " -f %s" % output_format if output_format else ""
        options += " -i %s" % per_thread_memory if per_thread_memory else ""

        options += " %s" % lastdb
        options += " %s" % query
        if output_format == "MAF":
            output_filename_list = self.split_filename(output)
            tab_filename = output + ".tab" if output_filename_list[-1] != "maf" else output_filename_list[0] + output_filename_list[1] + ".tab"
            options += " | tee %s | maf-convert tab > %s" % (output_format, tab_filename)
        else:
            options += " > %s" % output
        return options

    def lastal(self, lastdb, query, output, verbose=True,
               keep_preliminary_masking=True, mask_simple_repeats=True,
               output_format="MAF", per_thread_memory="4G"):

        options = self.parse_lastal_options(lastdb, query, output, verbose=verbose,
                                            keep_preliminary_masking=keep_preliminary_masking,
                                            mask_simple_repeats=mask_simple_repeats,
                                            output_format=output_format,
                                            per_thread_memory=per_thread_memory)

        self.execute(options=options, cmd="lastal")

    def convert_maf(self, maf_file, output_file, format="TAB"):

        options = " %s" % format
        options += " %s" % maf_file
        options += " > %s" % output_file

        self.execute(options=options, cmd="maf-convert")
