#!/usr/bin/env python
import os

from Routines.File import split_filename

from Tools.Abstract import Tool


class RAxML(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "raxml", path=path, max_threads=max_threads)

    def construct_tree_with_bootstrap(self,
                                      alignment_file,
                                      partition_file,
                                      output_prefix,
                                      model="GTRGAMMAI",
                                      bootstrap_replicates=100,
                                      outgroup=None,
                                      bootstrap_random_seed=12345,
                                      parsimony_random_seed=12345):
        options = " -f a"
        options += " -T %i" % self.threads
        options += " -x %s" % alignment_file
        options += " -n %s" % output_prefix
        options += " -o %s" % (outgroup if isinstance(outgroup, str) else ",".join(outgroup)) if outgroup else ""
        options += " -q %s" % partition_file
        options += " -# %s" % bootstrap_replicates
        options += " -x %i" % bootstrap_random_seed
        options += " -p %i" % parsimony_random_seed
        options += " -m %s" % model

        self.execute(options)
