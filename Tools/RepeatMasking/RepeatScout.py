#!/usr/bin/env python
import os

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import IdSet


class RepeatScout(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "RepeatScout", path=path, max_threads=max_threads)

    def predict_repeats(self, sequence, output, freq, l):
        """
        Usage:
        RepeatScout -sequence <seq> -output <out> -freq <freq> -l <l> [opts]
             -L # size of region to extend left or right (10000)
             -match # reward for a match (+1)
             -mismatch # penalty for a mismatch (-1)
             -gap  # penalty for a gap (-5)
             -maxgap # maximum number of gaps allowed (5)
             -maxoccurrences # cap on the number of sequences to align (10,000)
             -maxrepeats # stop work after reporting this number of repeats (10000)
             -cappenalty # cap on penalty for exiting alignment of a sequence (-20)
             -tandemdist # of bases that must intervene between two l-mers for both to be counted (500)
             -minthresh # stop if fewer than this number of l-mers are found in the seeding phase (3)
             -minimprovement # amount that a the alignment needs to improve each step to be considered progress (3)
             -stopafter # stop the alignment after this number of no-progress columns (100)
             -goodlength # minimum required length for a sequence to be reported (50)
             -maxentropy # entropy (complexity) threshold for an l-mer to be considered (-.7)
             -v[v[v[v]]] How verbose do you want it to be?  -vvvv is super-verbose.
        """

        options = " -sequence %s" % sequence
        options += " -output %s" % output
        options += " -freq %s" % str(freq)
        options += " -l %s" % str(l)

        self.cmd(options)