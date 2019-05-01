#!/usr/bin/env python
import os
#import shutil
#from collections import OrderedDict
#from RouToolPa.Tools.Filter import Cookiecutter, Trimmomatic, FaCut
from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Routines.Matplotlib import MatplotlibRoutines
from RouToolPa.Routines.Fastq import FastQRoutines
from RouToolPa.Routines.Annotations import AnnotationsRoutines
#from Parsers.FaCut import FaCutReport
#from Parsers.Coockiecutter import CoockiecutterReport
#from Parsers.Trimmomatic import TrimmomaticReport

#from RouToolPa.Collections.General import TwoLvlDict


class Pipeline(Tool, MatplotlibRoutines, FastQRoutines, AnnotationsRoutines):
    def __init__(self, max_threads=1, max_memory=10, workdir="./"):
        Tool.__init__(self, cmd="", max_threads=max_threads, max_memory=max_memory)
        self.workdir = workdir
        self.dirs = None

    @staticmethod
    def get_sample_list(sample_dir, sample_list=None):
        samples = []
        if sample_list:
            return [sample_list] if isinstance(sample_list, str) else sample_list
        else:
            dir_list = os.listdir(sample_dir)
            for directory in dir_list:
                if os.path.isdir("%s/%s" % (sample_dir, directory)):
                    samples.append(directory)

            return samples

    def init_dirs(self, dirs=None, workdir=None):

        working_directory = workdir if workdir else self.workdir if self.workdir else "./"

        self.recursive_mkdir(dir_dict=dirs if dirs else self.dirs,
                             out_dir=working_directory)
