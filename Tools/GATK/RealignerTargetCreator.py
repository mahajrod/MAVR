#!/usr/bin/env python
__author__ = 'mahajrod'

from Tools.GATK.Abstract import GATKTool

class RealignerTargetCreator(GATKTool):

    def find(self, reference, alignment, output="forIndelRealigner.intervals", known_indels_vcf=None,
             max_interval_size=None, min_reads_cov=None, mismatch_fraction=None, window_size=None):

        options = ""
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -o %s" % output
        options += " --known %s" % known_indels_vcf if known_indels_vcf else ""
        options += " -maxInterval %i" % max_interval_size if max_interval_size else ""
        options += " -minReads %i" % min_reads_cov if min_reads_cov else ""
        options += " -mismatch %i" % mismatch_fraction if mismatch_fraction else ""
        options += " -window %i" % window_size if window_size else ""

        self.execute(options, cmd="-T RealignerTargetCreator")