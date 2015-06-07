#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class ToolTemplate(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "tool_command", path=path, max_threads=max_threads)

    def tool_method(self, tool_options):
        options = ""
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        self.execute(options, cmd="custom command if necessary ex: tool_command subcommand")


if __name__ == "__main__":
    pass