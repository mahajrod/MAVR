__author__ = 'mahajrod'
import os
import re
import sys
import pickle

from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from CustomCollections.GeneralCollections import TwoLvlDict, SynDict, IdList, IdSet
from Routines.Functions import output_dict


class FastQRoutines():

    def __init__(self):
        pass

    @staticmethod
    def reverse_complement(in_file, out_file):
        with open(in_file, "r") as in_fd:
            with open(out_file, "w") as out_fd:
                for line in in_fd:
                    out_fd.write(line)
                    out_fd.write(str(Seq(in_fd.next().strip()).reverse_complement()))
                    out_fd.write("\n")
                    out_fd.write(in_fd.next())
                    out_fd.write(in_fd.next().strip()[::-1])
                    out_fd.write("\n")
