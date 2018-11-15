#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class BionanoSolve(Tool):
    """
    Requirements:
    1. python 2.7.5
    2. perl 5.14.X or 5.16.X
    3. R 3.1.2 or greater with R libraries data.table, igraph, intervals, MASS, parallel, XML
    (which requires a Linux system library libxml2), argparser
    4. glibc and gcc libraries (for older CPUs without AVX)
    5. minimum RAM size is 256 GB on at least 1 node, 32 GB on all nodes
    Running the assembly pipeline on a cluster requires a batch system that supports the Distributed
    Resource Management Application API (DRMAA). The submission hosts for the cluster also require
    the Python DRMAA module.
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "", path=path, max_threads=max_threads)

    def two_enzyme_scaffolding(self, parameter_file, enzyme1_map, enzyme2_map,
                               enzyme1_name, enzyme2_name, assembly_fasta,
                               ref_aligner_path, output_prefix, output_dir="./",
                               enzyme1_manual_cut=None, enzyme2_manual_cut=None):
        options = " "
        options += " -b1 %s" % enzyme1_map
        options += " -b2 %s" % enzyme2_map
        options += " -e1 %s" % enzyme1_name
        options += " -e2 %s" % enzyme2_name
        options += " -m1 %s" % enzyme1_manual_cut
        options += " -m2 %s" % enzyme2_manual_cut

        options += " -N %s" % assembly_fasta
        options += " -R %s" % ref_aligner_path
        options += " -O %s" % output_dir
        options += " -t %s/%s.tar" % (output_dir, output_prefix)
        options += " -s %s/%s.status" % (output_dir, output_prefix)

        options += " %s" % parameter_file

        self.execute(options, cmd="runTGH.R", intepreter="Rscript")

    def error_correction(self):
        pass
