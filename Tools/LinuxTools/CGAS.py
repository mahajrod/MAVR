#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from collections import Iterable
from subprocess import PIPE, Popen

from Tools.Abstract import Tool


class CGAS(Tool):
    def __init__(self, path="", max_threads=4):

        Tool.__init__(self, "", path=path, max_threads=max_threads)

    @staticmethod
    def _check_input(input_files):

        if not isinstance(input_files, Iterable):
            raise ValueError("%s is not iterable" % input_files)
        return input_files if isinstance(input_files, str) else " ".join(input_files)

    def cat(self, input_files, output=None, capture_output=False):

        options = " %s" % self._check_input(input_files)
        options += " > %s" % output if output else ""

        self.execute(options, cmd="cat", capture_output=capture_output)

    @staticmethod
    def _parse_grep_options(pattern, count_mode=False, use_regexp=False, num_lines_after_match=None,
                            num_lines_before_match=None, invert_match=False, ignore_case=False, whole_word_match=False,
                            max_matches=None, show_line_numbers=False, recursive_search=False):

        options = " -P" if use_regexp else ""
        options += " -c" if count_mode else ""
        options += " -A %i" % num_lines_after_match if num_lines_after_match else ""
        options += " -B %i" % num_lines_before_match if num_lines_before_match else ""
        options += " -v" if invert_match else ""
        options += " -i" if ignore_case else ""
        options += " -w" if whole_word_match else ""
        options += " -m %i" % max_matches if max_matches else ""
        options += " -n" if show_line_numbers else ""
        options += " -r" if recursive_search else ""
        options += " %s" % pattern

        return options

    def grep(self, pattern, input_files, output=None, count_mode=False, use_regexp=False, num_lines_after_match=None,
             num_lines_before_match=None, invert_match=False, ignore_case=False, whole_word_match=False,
             max_matches=None, show_line_numbers=False, recursive_search=False, capture_output=False):

        options = self._parse_grep_options(pattern, count_mode=count_mode, use_regexp=use_regexp,
                                           num_lines_after_match=num_lines_after_match,
                                           num_lines_before_match=num_lines_before_match, invert_match=invert_match,
                                           ignore_case=ignore_case, whole_word_match=whole_word_match,
                                           max_matches=max_matches, show_line_numbers=show_line_numbers,
                                           recursive_search=recursive_search)
        options += " %s" % self._check_input(input_files)
        options += " > %s" % output if output else ""

        self.execute(options, cmd="grep", capture_output=capture_output)

    @staticmethod
    def _parse_awk_options(awk_code, separator=None):

        options = " -F'%s'" % separator
        options += " '%s'" % awk_code

        return options

    def awk(self, awk_code, input_files, output=None, separator=None, capture_output=False):

        options = self._parse_awk_options(awk_code, separator=separator)
        options += " %s" % self._check_input(input_files)
        options += " > %s" % output if output else ""

        self.execute(options, cmd="awk", capture_output=capture_output)

    @staticmethod
    def _parse_sed_options(sed_string, use_regexp=False, disable_print_of_each_line=False):

        options = " -r" if use_regexp else ""
        options += " -n" if disable_print_of_each_line else ""
        options += " '%s'" % sed_string

        return options

    def sed(self, sed_string, input_files, output=None, use_regexp=False, disable_print_of_each_line=False,
            capture_output=False):

        options = self._parse_sed_options(sed_string, use_regexp=use_regexp,
                                          disable_print_of_each_line=disable_print_of_each_line)

        options += " %s" % self._check_input(input_files)
        options += " > %s" % output if output else ""

        self.execute(options, cmd="sed", capture_output=capture_output)

    def cgas(self, input_files, grep_pattern=None, awk_code=None, sed_string=None, output=None, count_mode=False,
             grep_use_regexp=False, num_lines_after_match=None, num_lines_before_match=None, invert_match=False,
             ignore_case=False, whole_word_match=False, max_matches=None, show_line_numbers=False,
             recursive_search=False, separator=None, sed_use_regexp=False, disable_print_of_each_line=False,
             capture_output=False):

        options = self._check_input(input_files)

        if grep_pattern:
            options += " | grep %s" % self._parse_grep_options(grep_pattern, count_mode=count_mode,
                                                               use_regexp=grep_use_regexp,
                                                               num_lines_after_match=num_lines_after_match,
                                                               num_lines_before_match=num_lines_before_match,
                                                               invert_match=invert_match,
                                                               ignore_case=ignore_case,
                                                               whole_word_match=whole_word_match,
                                                               max_matches=max_matches,
                                                               show_line_numbers=show_line_numbers,
                                                               recursive_search=recursive_search)
        if awk_code:
            options += " | awk %s" % self._parse_awk_options(awk_code, separator=separator)

        if sed_string:
            options += " | sed %s" % self._parse_sed_options(sed_string, use_regexp=sed_use_regexp,
                                                             disable_print_of_each_line=disable_print_of_each_line)

        options += " "

        options += " > %s" % output if output else ""
        self.execute(options, cmd="cat", capture_output=capture_output)


if __name__ == "__main__":
    pass