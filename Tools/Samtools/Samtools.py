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
        self.bam_flags = {
                          "read_paired": 1,
                          "read_mapped_in_proper_pair": 2,
                          "read_unmapped": 4,
                          "mate_unmapped": 8,
                          "read_reverse_strand": 16,
                          "mate_reverse_strand": 32,
                          "first_in_pair": 64,
                          "second_in_pair": 128,
                          "not_primary_alignment": 256,
                          "read_fails_platform/vendor_quality_checks": 512,
                          "read_is_PCR_or_optical_duplicate": 1024,
                          "supplementary_alignment": 2048
                          }

    def rmdup(self, input_bam, output_bam, remove_dup_for_se_reads=False, treat_both_pe_and_se_reads=False):

        options = " -s" if remove_dup_for_se_reads else ""
        options += " -S" if treat_both_pe_and_se_reads else ""
        options += " %s" % input_bam
        options += " %s" % output_bam

        self.execute(options, cmd="samtools rmdup")

    def parse_view_options(self, include_header_in_output=True, output_uncompressed_bam=False,
                           output_bam=False, white_flag_value=None, black_flag_value=None,
                           bed_file_with_regions_to_output=None):
        options = " -@ %i" % self.threads
        options += " -h" if include_header_in_output else ""
        options += " -u" if output_uncompressed_bam else ""
        options += " -b" if output_bam else ""
        options += " -f %i" % white_flag_value if white_flag_value else ""
        options += " -F %i" % black_flag_value if black_flag_value else ""
        options += " -L %s" % bed_file_with_regions_to_output if bed_file_with_regions_to_output else ""
        return options

    def view(self, input_file, output_file, include_header_in_output=True, output_uncompressed_bam=False,
             output_bam=False, white_flag_value=None, black_flag_value=None, bed_file_with_regions_to_output=None):

        options = self.parse_view_options(include_header_in_output=include_header_in_output,
                                          output_uncompressed_bam=output_uncompressed_bam,
                                          output_bam=output_bam, white_flag_value=white_flag_value,
                                          black_flag_value=black_flag_value,
                                          bed_file_with_regions_to_output=bed_file_with_regions_to_output)

        options += " -o %s" % output_file
        options += " %s" % input_file

        self.execute(options, "samtools view")

    def parse_sort_options(self, temp_file_prefix="temp_bam", sort_by_name=False, max_memory_per_thread="1G"):
        options = " -@ %i" % self.threads
        options += " -n" if sort_by_name else ""
        options += " -m %s" % max_memory_per_thread if max_memory_per_thread else ""
        options += " -T %s" % temp_file_prefix
        return options

    def sort(self, input_bam, output_bam, temp_file_prefix="temp_bam", sort_by_name=False, max_memory_per_thread="1G"):
        options = self.parse_sort_options(temp_file_prefix=temp_file_prefix, sort_by_name=sort_by_name,
                                          max_memory_per_thread=max_memory_per_thread)
        options += " -o %s" % output_bam
        options += " %s" % input_bam

        self.execute(options, "samtools sort")

    def get_info_from_flags(self, flag_list):
        info_field = 0
        for flag in flag_list:
            info_field |= self.bam_flags[flag]

        return info_field

    def prepare_bam_for_read_extraction(self, input_bam, output_bam, temp_file_prefix="temp_bam",
                                        max_memory_per_thread="1G"):
        black_list_flags = ["not_primary_alignment", "supplementary_alignment"]
        info_field = self.get_info_from_flags(black_list_flags)

        view_options = " -b -u -F %i %s" % (info_field, input_bam)
        sort_options = " -n -m %s -T %s -o %s -" % (max_memory_per_thread, temp_file_prefix, output_bam)
        cmd = "samtools view %s | samtools sort %s" % (view_options, sort_options)
        self.execute(options="", cmd=cmd)

    def faidx(self, fasta_file):

        options = " %s" % fasta_file

        self.execute(options, cmd="samtools faidx")

    def index(self, bam_file):

        options = " %s" % bam_file

        self.execute(options, cmd="samtools index")


class SamtoolsV0(SamtoolsV1, Tool):
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

    def sort(self, input_bam, output_prefix, max_memory=3000000000):

        options = " %s" % input_bam
        options += " %s" % output_prefix
        options += " -m %i" % max_memory

        self.execute(options, "samtools sort")

