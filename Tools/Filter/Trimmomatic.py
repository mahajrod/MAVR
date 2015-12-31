__author__ = 'mahajrod'
import os
import pyparsing as pyp
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

        options += " MINLEN:%i" % min_length if min_length else ""

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

        options += " > %s 2>&1" % logfile if logfile else ""

        self.execute(options=options)

    @staticmethod
    def parse_log(log_file):

        input = pyp.Literal("Input Read Pairs: ")
        input_number = pyp.Word(pyp.nums)
        both = pyp.Literal("Both Surviving: ")
        both_number = pyp.Word(pyp.nums)

        str1 = pyp.Literal("(")
        both_percent = pyp.Word(pyp.nums + ".-")
        forward = pyp.Literal("%) Forward Only Surviving: ")
        forward_number = pyp.Word(pyp.nums)
        forward_percent = pyp.Word(pyp.nums + ".-")
        reverse = pyp.Literal("%) Reverse Only Surviving: ")
        reverse_number = pyp.Word(pyp.nums)
        reverse_percent = pyp.Word(pyp.nums + ".-")
        dropped = pyp.Literal("%) Dropped: ")
        dropped_number = pyp.Word(pyp.nums)
        dropped_percent = pyp.Word(pyp.nums + ".-")
        end = pyp.Literal("%)")

        for token in both_number, forward_number, reverse_number, dropped_number, input_number:
            token.setParseAction(lambda s, l, t: [int(t[0])])
        for token in both_percent, forward_percent, reverse_percent, dropped_percent:
            token.setParseAction(lambda s, l, t: [float(t[0])])

        both_number.setParseAction(lambda s, l, t: [int(t[0])])

        pattern = input + input_number + both + both_number + str1 + both_percent + forward + forward_number + str1 + \
                  forward_percent + reverse + reverse_number + str1 + reverse_percent + dropped + dropped_number + \
                  str1 + dropped_percent + end

        with open(log_file, "r") as log_fd:
            for line in log_fd:
                if line[:17] == "Input Read Pairs:":
                    parsed_list = pattern.parseString(line.strip())

        return {
                "input":                parsed_list[1],
                "both_surviving":       parsed_list[3],
                "both_surviving,%":     parsed_list[5],
                "forward_only_surviving":    parsed_list[7],
                "forward_only_surviving,%":  parsed_list[9],
                "reverse_only_surviving":    parsed_list[11],
                "reverse_only_surviving,%":  parsed_list[13],
                "dropped_only_surviving":    parsed_list[15],
                "dropped_only_surviving,%":  parsed_list[17]
                }



