#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from collections import Iterable
from subprocess import PIPE, Popen

from Tools.Abstract import Tool


class Valgrind(Tool):
    """
    Memory usage tester and profiler
    """
    def __init__(self, path="", max_threads=4):

        Tool.__init__(self, "valgrind", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(cmd, cmd_log=None, massif_mode=False,  count_all_memory=False,
                      valgrind_log_file=None, massif_report_file=None):
        options = ""
        options += " --tool=massif" if massif_mode else ""
        options += " --pages-as-heap=yes" if count_all_memory else ""
        options += " --log-file=%s" % valgrind_log_file if valgrind_log_file else ""
        options += " --massif-out-file=%s" % massif_report_file if massif_report_file else ""
        options += " -b %s" % cmd
        options += " > %s 2>&1" % cmd_log if cmd_log else ""

        return options

    def ms_print(self, massif_report, ms_print_report):

        options = " %s" % massif_report
        options += " > %s" % ms_print_report

        self.execute(options=options, cmd="ms_print")

    def measure_ram_usage(self, cmd, output_prefix, count_all_memory=True):

        cmd_log = "%s.cmd.log" % output_prefix
        valgrind_log = "%s.valgrind.log" % output_prefix
        massif_report = "%s.valgrind.massif.out" % output_prefix
        ms_print_report = "%s.valgrind.massif.ms_print.out" % output_prefix

        options = self.parse_options(cmd, cmd_log=cmd_log, massif_mode=True,  count_all_memory=count_all_memory,
                                     valgrind_log_file=valgrind_log, massif_report_file=massif_report)

        self.execute(options=options)
        self.ms_print(massif_report, ms_print_report)

if __name__ == "__main__":
    pass
