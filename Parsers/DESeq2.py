#!/usr/bin/python2
import os
import re
from collections import OrderedDict
from math import sqrt

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from Parsers.Abstract import Metadata, Header


class RecordPWC:
    def __init__(self, gene, basemean, log2foldchange, lfcSE, stat, pvalue, padj,):
        self.gene = gene                              # str
        self.basemean = basemean                      # float
        self.log2foldchange = log2foldchange          # float
        self.lfcSE = lfcSE                            # float
        self.stat = stat                              # float
        self.pvalue = pvalue                          # float
        self.padj = padj                              # float


class MetadataPWC(OrderedDict, Metadata):
    def __init__(self):
        OrderedDict.__init__(self)
        Metadata.__init__(self)


class HeaderPWC(list, Header):

    def __str__(self):
        return "\t".join(self)


class CollectionPWC(OrderedDict):
    def __init__(self, metadata=None, record_dict=None, header=None, pwc_file=None, from_file=True):
        OrderedDict.__init__(self)
        if from_file:
            self.metadata = Metadata()
            with open(pwc_file, "r") as fd:
                self.header = HeaderPWC(fd.readline().strip().split())
                for line in fd:
                    record = self.add_record(line)
                    self[record.gene] = record
        else:
            self.metadata = metadata
            self.header = header
            for gene in record_dict:
                self[gene] = record_dict[gene]

    @staticmethod
    def add_record(line):
        line_list = line.strip().split("\t")
        return RecordPWC(line_list[0].strip("\""), float(line_list[1]), float(line_list[2]), float(line_list[3]),
                         float(line_list[4]), float(line_list[5]), float(line_list[6]))

