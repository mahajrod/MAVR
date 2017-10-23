#!/usr/bin/env python
__author__ = 'mahajrod'

import pyparsing as pyp


class TrimmomaticReport:
    def __init__(self, trimmomatic_report, input_is_se=False):

        if input_is_se:
            input = pyp.Literal("Input Reads: ")
            input_number = pyp.Word(pyp.nums)
            surv = pyp.Literal("Surviving: ")
            surv_number = pyp.Word(pyp.nums)

            str1 = pyp.Literal("(")
            surv_percent = pyp.Word(pyp.nums + ".-")
            dropped = pyp.Literal("%) Dropped: ")
            dropped_number = pyp.Word(pyp.nums)
            dropped_percent = pyp.Word(pyp.nums + ".-")
            end = pyp.Literal("%)")

            for token in surv_number, dropped_number, input_number:
                token.setParseAction(lambda s, l, t: [int(t[0])])
            for token in surv_percent, dropped_percent:
                token.setParseAction(lambda s, l, t: [float(t[0])])

            pattern = input + input_number + surv + surv_number + str1 + surv_percent + dropped + dropped_number + \
                      str1 + dropped_percent + end

            with open(trimmomatic_report, "r") as log_fd:
                for line in log_fd:
                    if line[:12] == "Input Reads:":
                        parsed_list = pattern.parseString(line.strip())

            self.stats = {
                          "input":                   parsed_list[1],
                          "surviving":               parsed_list[3],
                          "surviving,%":             parsed_list[5],
                          "dropped":                 parsed_list[7],
                          "dropped,%":               parsed_list[9]
                         }

        else:

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

            with open(trimmomatic_report, "r") as log_fd:
                for line in log_fd:
                    if line[:17] == "Input Read Pairs:":
                        parsed_list = pattern.parseString(line.strip())

            self.stats = {
                          "input":                     parsed_list[1],
                          "both_surviving":            parsed_list[3],
                          "both_surviving,%":          parsed_list[5],
                          "forward_only_surviving":    parsed_list[7],
                          "forward_only_surviving,%":  parsed_list[9],
                          "reverse_only_surviving":    parsed_list[11],
                          "reverse_only_surviving,%":  parsed_list[13],
                          "dropped":                   parsed_list[15],
                          "dropped,%":                 parsed_list[17]
                         }


