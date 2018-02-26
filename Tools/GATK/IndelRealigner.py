#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.GATK.Abstract import GATKTool


class IndelRealigner(GATKTool):

    # "find" method was renamed to realign
    def realign(self, reference, alignment, corrected_alignment, target_intervals="forIndelRealigner.intervals",
                known_indels_vcf=None, model=None, lod_threshold=None, entropy_threshold=None, max_cons=None,
                max_size_for_movement=None, max_pos_move=None, max_reads_for_cons=None,
                max_reads_for_realignment=None, max_reads_in_memory=None, no_original_tags=False,
                nway_out=False, default_base_qualities=None):

        options = " -nt %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " -o %s" % corrected_alignment
        options += " -targetIntervals %s" % target_intervals
        options += " --known %s" % known_indels_vcf if known_indels_vcf else ""
        options += " -model %s" % model if model else ""
        options += " -LOD %i" % lod_threshold if lod_threshold else ""
        options += " -entropy %i" % entropy_threshold if entropy_threshold else ""
        options += " --maxConsensuses %i" % max_cons if max_cons else ""
        options += " -maxisize %i" % max_size_for_movement if max_size_for_movement else ""
        options += " -maxPosMove %i" % max_pos_move if max_pos_move else ""
        options += " -greedy %i" % max_reads_for_cons if max_reads_for_cons else ""
        options += " -maxReads %i" % max_reads_for_realignment if max_reads_for_realignment else ""
        options += " -maxinMemory %i" % max_reads_in_memory if max_reads_in_memory else ""

        options += " --nWayOut" if nway_out else ""
        options += " -noTags" if no_original_tags else ""
        options += " --defaultBaseQualities %i" % default_base_qualities if default_base_qualities else ""

        self.execute(options, cmd="-T IndelRealigner")

if __name__ == "__main__":
    import os
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc"
    os.chdir(workdir)
    #gatk_dir =
    reference = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    alignment = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/001/alignment_tmap_Alsu24mc/001_trimmed_sorted_rm_pcr_chrom.bam"
    corrected_alignment = "001_trimmed_sorted_rm_pcr_chrom_corrected.bam"
    IndelRealigner = IndelRealigner(jar_path="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0")
    IndelRealigner.find(reference, alignment, corrected_alignment)
