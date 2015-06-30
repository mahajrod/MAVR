#!/usr/bin/env python

from Tools.Abstract import Tool


class SamtoolsV1(Tool):
    """
    Class for samtools 1.0+
    Several subcommands are not implemented
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "samtools", path=path, max_threads=max_threads)

        # bam/sam flag values:
        self.bam_flag_values = {"unaligned": 4,
                                "supplementary_alignment": 2048,
                                "non_primary_alignment": 256
                                }

    def rmdup(self, input_bam, output_bam, remove_dup_for_se_reads=False, treat_both_pe_and_se_reads=False):

        options = " -s" if remove_dup_for_se_reads else ""
        options += " -S" if treat_both_pe_and_se_reads else ""
        options += " %s" % input_bam
        options += " %s" % output_bam

        self.execute(options, cmd="samtools rmdup")

    def view(self, input_file, output_file, include_header_in_output=True, output_uncompressed_bam=False,
             output_bam=False, white_flag_value=None, black_flag_value=None, bed_file_with_regions_to_output=None):

        options = " -@ %i" % self.threads
        options += " -h" if include_header_in_output else ""
        options += " -u" if output_uncompressed_bam else ""
        options += " -b" if output_bam else ""
        options += " -o %s" % output_file
        options += " -f %i" % white_flag_value if white_flag_value else ""
        options += " -F %i" % black_flag_value if black_flag_value else ""
        options += " -L %s" % bed_file_with_regions_to_output if bed_file_with_regions_to_output else ""
        options += " %s" % input_file

        self.execute(options, "samtools view")

    def sort(self, input_bam, output_bam, temp_file_prefix="temp_bam"):

        options = " -@ %i" % self.threads
        options += " -o %s" % output_bam
        options += " -T %s" % temp_file_prefix
        options += " %s" % input_bam
        self.execute(options, "samtools sort")

    def faidx(self, fasta_file):

        options = " %s" % fasta_file

        self.execute(options, cmd="samtools faidx")

    def index(self, bam_file):

        options = " %s" % bam_file

        self.execute(options, cmd="samtools index")


class SamtoolsV0(Tool, SamtoolsV1):
    """
    Class for samtools 0.1.19+
    Several subcommands are not implemented
    """

    def view(self, input_file, output_file, sam_input=False, include_header_in_output=True,
             output_uncompressed_bam=False, output_bam=False, white_flag_value=None, black_flag_value=None,
             bed_file_with_regions_to_output=None):

        options = " -h" if include_header_in_output else ""
        options += " -S" if sam_input else ""
        options += " -u" if output_uncompressed_bam else ""
        options += " -b" if output_bam else ""
        options += " -o %s" % output_file
        options += " -f %i" % white_flag_value if white_flag_value else ""
        options += " -F %i" % black_flag_value if black_flag_value else ""
        options += " -L %s" % bed_file_with_regions_to_output if bed_file_with_regions_to_output else ""
        options += " %s" % input_file

        self.execute(options, "samtools view")

    def sort(self, input_bam, output_prefix):

        options = " %s" % input_bam
        options += " %s" % output_prefix

        self.execute(options, "samtools sort")

