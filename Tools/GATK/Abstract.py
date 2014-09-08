#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import Tool
from Routines.Functions import check_path


class GATKTool(Tool):

    def __init__(self, path="", max_threads=4, jar_path="", jar="GenomeAnalysisTK.jar",
                 max_memory="1g"):
        self.cmd = "java"
        self.path = check_path(path)
        self.max_threads = max_threads
        self.jar_path = check_path(jar_path)
        self.jar = jar
        self.max_memory = max_memory