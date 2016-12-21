#!/usr/bin/env python
import os
import shutil

from collections import OrderedDict

from Tools.Filter import Cookiecutter, Trimmomatic, FaCut
from Routines import FileRoutines
from Parsers.FaCut import FaCutReport
from Parsers.Coockiecutter import CoockiecutterReport
from Parsers.Trimmomatic import TrimmomaticReport

from CustomCollections.GeneralCollections import TwoLvlDict


class Pipeline:
    def __init__(self):
        pass

    @staticmethod
    def get_sample_list(samples_directory):
        samples = sorted(os.listdir(samples_directory))
        sample_list = []
        for sample in samples:
            if os.path.isdir("%s/%s" % (samples_directory, sample)):
                sample_list.append(sample)
        return sample_list

    @staticmethod
    def combine_fastq_files(samples_directory, sample, output_directory):
        sample_dir = "%s/%s/" % (samples_directory, sample)
        filetypes, forward_files, reverse_files = FileRoutines.make_lists_forward_and_reverse_files(sample_dir)
        if len(filetypes) == 1:
            if "fq.gz" in filetypes:
                command = "zcat"
            elif "fq.bz2" in filetypes:
                command = "bzcat"
            else:
                command = "cat"

            os.system("%s %s > %s/%s_1.fq" % (command, " ".join(forward_files), output_directory, sample))
            os.system("%s %s > %s/%s_2.fq" % (command, " ".join(reverse_files), output_directory, sample))
        else:
            raise IOError("Extracting from mix of archives in not implemented yet")
