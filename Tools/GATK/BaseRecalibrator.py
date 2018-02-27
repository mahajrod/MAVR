#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool


class BaseRecalibrator(JavaTool):

    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T BaseRecalibrator", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def get_recalibration_table(self, reference, alignment, output_table, known_sites_vcf):

        # TODO: add rest of  options
        options = ""
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -nct %i" % self.threads
        options += " -knownSites %s" % known_sites_vcf if isinstance(known_sites_vcf, str) else " -knownSites ".join(known_sites_vcf)
        options += " -o %s" % output_table

        self.execute(options)

if __name__ == "__main__":
    import os
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc"
    os.chdir(workdir)
    #gatk_dir =
    reference = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    alignment = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc/001_trimmed_sorted_rm_pcr_chrom_corrected.bam"
    corrected_alignment = "001_trimmed_sorted_rm_pcr_chrom_corrected.bam"
    BaseRecalibrator = BaseRecalibrator(jar_path="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0")
    BaseRecalibrator.get_recalibration_table(reference, alignment)
