#/usr/bin/env python
import os
import shutil
from collections import OrderedDict
from RouToolPa.Tools.Alignment import Bowtie2, BWA, Novoalign
from RouToolPa.Tools.Samtools import SamtoolsV1
from RouToolPa.Tools.Picard import MarkDuplicates, AddOrReplaceReadGroups
from Pipelines.Abstract import Pipeline


class AlignmentPipeline(Pipeline):

    def __init__(self, max_threads=1, max_memory=10, BWA_dir="", BWA_binary=None, bowtie2_dir="", bowtie2_binary=None,
                 Picard_dir=""):
        Pipeline.__init__(self, max_threads=max_threads, max_memory=max_memory)

        self.BWA_dir = BWA_dir
        self.bowtie2_dir = bowtie2_dir
        self.Picard_dir = Picard_dir

        self.BWA_binary = BWA_binary
        self.bowtie2_binary = bowtie2_binary

    def prepare_dirs(self, sample_list, outdir="./", ):
        dir_dict = OrderedDict()

        for sample in sample_list:
            dir_dict[sample] = OrderedDict()
                
        self.recursive_mkdir(dir_dict, out_dir=outdir)

    def init_tools(self, threads=None):
        for tool in (MarkDuplicates,
                     AddOrReplaceReadGroups,
                     BWA,
                     Bowtie2,
                     Novoalign,
                     SamtoolsV1):

            tool.threads = threads if threads else self.threads
            tool.max_memory = "%ig" % self.max_memory

        BWA.path = self.BWA_dir
        Bowtie2.path = self.bowtie2_dir
        MarkDuplicates.jar_path = self.Picard_dir
        AddOrReplaceReadGroups.jar_path = self.Picard_dir

        if self.BWA_binary:
            BWA.cmd = BWA.cmd
        if self.bowtie2_binary:
            Bowtie2.cmd = self.bowtie2_binary

    def align(self, sample_dir, reference_index, aligner="bwa", sample_list=None, outdir="./",
              quality_score_type="phred33", read_suffix="", read_extension="fastq",
              alignment_format="bam", threads=None, mark_duplicates=True, platform="Illumina",
              add_read_groups_by_picard=False, gzipped_reads=False):

        self.init_tools(threads=threads)

        samples = self.get_sample_list(sample_dir, sample_list=sample_list)

        self.prepare_dirs(samples, outdir=outdir)

        if aligner == "bowtie2":
            aligner_tool = Bowtie2
        elif aligner == "bwa":
            aligner_tool = BWA
        else:
            raise ValueError("")

        for sample in samples:
            read_prefix = "%s/%s/%s%s" % (sample_dir, sample, sample, read_suffix)
            forward_reads = "%s_1.%s%s" % (read_prefix, read_extension, ".gz" if gzipped_reads else "")
            reverse_reads = "%s_2.%s%s" % (read_prefix, read_extension, ".gz" if gzipped_reads else "")

            output_prefix = "%s/%s/%s" % (outdir, sample, sample)

            raw_alignment = "%s.%s" % (output_prefix, alignment_format)
            final_alignment = "%s.mkdup.%s" % (output_prefix, alignment_format)

            duplicates_stat_file = "%s.duplicates.stat" % output_prefix
            coverage_file = "%s.coverage.bed" % output_prefix

            sorted_alignment_picard_groups = None

            aligner_tool.align(reference_index, forward_reads_list=forward_reads, reverse_reads_list=reverse_reads,
                               unpaired_reads_list=None, quality_score=quality_score_type, output_prefix=output_prefix,
                               output_format=alignment_format,
                               read_group_name=sample,
                               PU="x",
                               SM=sample,
                               platform=platform,
                               LB="x",
                               sort_by_coordinate=True,
                               sort_by_name=False,
                               max_per_sorting_thread_memory=str(max(int(self.max_memory/self.threads), 1)) + "G")

            if add_read_groups_by_picard:
                sorted_alignment_picard_groups = "%s.picard_groups.%s" % (output_prefix, alignment_format)
                AddOrReplaceReadGroups.add_read_groups(raw_alignment, sorted_alignment_picard_groups,
                                                       RGID=sample, RGLB=sample, RGPL=platform,
                                                       RGSM=sample, RGPU=sample)

            if alignment_format == "bam":
                SamtoolsV1.index(sorted_alignment_picard_groups if sorted_alignment_picard_groups else raw_alignment)

            if mark_duplicates:
                MarkDuplicates.run(sorted_alignment_picard_groups if sorted_alignment_picard_groups else raw_alignment,
                                   final_alignment,
                                   duplicates_stat_file)
                os.remove(raw_alignment)
                os.remove(raw_alignment + ".bai")
                if alignment_format == "bam":
                    SamtoolsV1.index(final_alignment)
