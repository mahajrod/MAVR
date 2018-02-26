#!/usr/bin/env python
__author__ = 'mahajrod'

from Tools.Abstract import Tool
from Tools.GATK.Abstract import GATKTool


class RealignerTargetCreator(GATKTool, Tool):

    def create(self, reference, alignment, output="forIndelRealigner.intervals", known_indels_vcf=None,
               max_interval_size=None, min_reads_cov=None, mismatch_fraction=None, window_size=None,
               default_base_qualities=None):

        options = " -nt %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -o %s" % output
        options += " --known %s" % known_indels_vcf if known_indels_vcf else ""
        options += " -maxInterval %i" % max_interval_size if max_interval_size else ""
        options += " -minReads %i" % min_reads_cov if min_reads_cov else ""
        options += " -mismatch %i" % mismatch_fraction if mismatch_fraction else ""
        options += " -window %i" % window_size if window_size else ""
        options += " --defaultBaseQualities %i" % default_base_qualities if default_base_qualities else ""

        self.execute(options, cmd="-T RealignerTargetCreator")

if __name__ == "__main__":
    import os
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc"
    os.chdir(workdir)
    #gatk_dir =
    reference = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    alignment = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc/001_trimmed_sorted_rm_pcr_chrom.bam"
    RealignerTargetCreator = RealignerTargetCreator(jar_path="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0")
    RealignerTargetCreator.create(reference, alignment)