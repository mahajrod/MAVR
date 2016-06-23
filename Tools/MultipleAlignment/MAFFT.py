#!/usr/bin/env python
import os

from Routines import FileRoutines

from Tools.Abstract import Tool


class MAFFT(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "mafft", path=path, max_threads=max_threads)

    def align(self,
              sequence_file, output, gap_open_penalty=None, offset=None, maxiterate=None,
              quiet=False, mode="globalpair", anysymbol=False):
        # TODO: add rest of options

        options = " --thread %i" % self.threads
        options += " --op %f" % gap_open_penalty if gap_open_penalty is not None else ""
        options += " --ep %f" % offset if offset is not None else ""
        options += " --maxiterate %i" % maxiterate if maxiterate is not None else ""
        options += " --quiet" if quiet else ""
        options += " --%s" % mode
        options += " --anysymbol" if anysymbol else ""

        options += " %s" % sequence_file
        options += " > %s" % output

        self.execute(options)

    def parallel_align(self, list_of_files, output_directory, output_suffix="alignment", gap_open_penalty=None,
                       offset=None, maxiterate=None, quiet=False, mode="globalpair", number_of_processes=1,
                       anysymbol=False):
        # TODO: add rest of options

        options = " --thread %i" % self.threads
        options += " --op %f" % gap_open_penalty if gap_open_penalty is not None else ""
        options += " --ep %f" % offset if offset is not None else ""
        options += " --maxiterate %i" % maxiterate if maxiterate is not None else ""
        options += " --quiet" if quiet else ""
        options += " --%s" % mode
        options += " --anysymbol" if anysymbol else ""
        options_list = []
        for filename in list_of_files:
            basename = FileRoutines.split_filename(filename)[1]
            op = options
            op += " %s" % filename
            op += " > %s/%s.fasta" % (output_directory,
                                      ("%s_%s" % (basename, output_suffix)) if output_suffix else basename)
            options_list.append(op)

        self.parallel_execute(options_list, threads=number_of_processes)