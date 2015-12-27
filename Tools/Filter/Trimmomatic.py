__author__ = 'mahajrod'
import os

from Tools.Abstract import JavaTool


class Trimmomatic(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", jar=None):
        jar = "trimmomatic-0.35.jar" if jar is None else jar
        JavaTool.__init__(self, jar, java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=None, timelog="trimmomatic.time.log")

    def parse_options(self, left_reads, output_prefix, output_extension="fq", right_reads=None, adapters_file=None,
                      mismatch_number=2, pe_reads_score=30, se_read_score=10, min_adapter_len=1,
                      sliding_window_size=None, average_quality_threshold=15,
                      leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
                      crop_length=None, head_crop_length=None, min_length=50, base_quality="phred33"):

        options = " PE" if right_reads else " SE"
        options += " -threads %i" % self.threads
        options += " -%s" % base_quality
        options += " %s" % left_reads
        options += " %s" % right_reads if right_reads else ""

        options += " %s_1.pe.%s" % (output_prefix, output_extension)
        options += " %s_1.se.%s" % (output_prefix, output_extension)
        options += " %s_2.pe.%s" % (output_prefix, output_extension)
        options += " %s_2.se.%s" % (output_prefix, output_extension)

        options += " ILLUMINACLIP:%s:%i:%i:%i:%i" % (adapters_file, mismatch_number, pe_reads_score,
                                                     se_read_score, min_adapter_len) if adapters_file else ""

        options += " SLIDINGWINDOW:%s:%s" % (sliding_window_size, average_quality_threshold) \
            if sliding_window_size else ""

        options += " LEADING:%i" % leading_base_quality_threshold if leading_base_quality_threshold else ""
        options += " TRAILING:%i" % trailing_base_quality_threshold if trailing_base_quality_threshold else ""

        options += " CROP:%i" % crop_length if crop_length else ""
        options += " HEADCROP:%i" % head_crop_length if head_crop_length else ""

        options += "MINLEN:%i" % min_length if min_length else ""

        return options

    def filter(self, left_reads, output_prefix, output_extension="fq", right_reads=None, adapters_file=None,
               mismatch_number=2, pe_reads_score=30, se_read_score=10, min_adapter_len=1,
               sliding_window_size=None, average_quality_threshold=15,
               leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
               crop_length=None, head_crop_length=None, min_length=50, logfile="trimmomatic.log",
               base_quality="phred33"):

        options = self.parse_options(left_reads, output_prefix, output_extension=output_extension,
                                     right_reads=right_reads, adapters_file=adapters_file,
                                     mismatch_number=mismatch_number, pe_reads_score=pe_reads_score,
                                     se_read_score=se_read_score, min_adapter_len=min_adapter_len,
                                     sliding_window_size=sliding_window_size, average_quality_threshold=average_quality_threshold,
                                     leading_base_quality_threshold=leading_base_quality_threshold,
                                     trailing_base_quality_threshold=trailing_base_quality_threshold,
                                     crop_length=crop_length, head_crop_length=head_crop_length, min_length=min_length,
                                     base_quality=base_quality)

        options += " > %s" % logfile if logfile else ""

        self.execute(options=options)
