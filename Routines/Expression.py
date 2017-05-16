#!/usr/bin/env python
__author__ = 'mahajrod'
import numpy as np
from CustomCollections.GeneralCollections import IdList


class ExpressionRoutines:
    def __init__(self):
        pass

    @staticmethod
    def divide_counts_by_base_level(input_file, output_prefix, base_level,
                                    separator="\t", verbose=True, secondary_base_lvl=None):
        output_file = "%s.divided_by_baselvl" % output_prefix
        zero_base_lvl_list = IdList()
        zero_both_base_lvls_list = IdList()
        zero_base_lvl_list_file = "%s.zero_base_lvl.ids" % output_prefix
        zero_both_base_lvls_list_file = "%s.zero_base_lvls.ids" % output_prefix
        with open(input_file, "r") as in_fd:
            header = in_fd.readline()
            header_list = header.strip().split(separator)
            data_base_level_index = header_list.index(base_level) - 1
            if secondary_base_lvl:
                data_secondary_base_level_index = header_list.index(secondary_base_lvl) - 1
            with open(output_file, "w") as out_fd:
                out_fd.write(header)
                for line in in_fd:
                    tmp_line = line.strip().split(separator)
                    data = np.array(map(float, tmp_line[1:]))
                    if data[data_base_level_index] == 0:
                        zero_base_lvl_list.append(tmp_line[0])
                        if not secondary_base_lvl:
                            if verbose:
                                print("Zero base level(%s) for %s...Skipping..." % (base_level, tmp_line[0]))
                            continue
                    if secondary_base_lvl:
                        if data[data_secondary_base_level_index] == 0:
                            zero_both_base_lvls_list.append(tmp_line[0])
                            if verbose:
                                print("Both base levels are zero (%s, %s) for %s...Skipping..." % (base_level,
                                                                                                   secondary_base_lvl,
                                                                                                   tmp_line[0]))
                            continue

                        data /= data[data_base_level_index] if data[data_base_level_index] != 0 else data[data_secondary_base_level_index]
                    else:
                        data /= data[data_base_level_index]
                    output_string = tmp_line[0] + "\t"
                    output_string += "\t".join(map(str, data))
                    output_string += "\n"
                    out_fd.write(output_string)

        zero_base_lvl_list.write(zero_base_lvl_list_file)
        zero_both_base_lvls_list.write(zero_both_base_lvls_list_file)

    @staticmethod
    def divide_counts_by_max_level(input_file, output_prefix,
                                   separator="\t", verbose=True,):
        output_file = "%s.divided_by_maxlvl" % output_prefix
        zero_max_lvl_list = IdList()

        zero_max_lvl_list_file = "%s.zero_max_lvl.ids" % output_prefix

        with open(input_file, "r") as in_fd:
            header = in_fd.readline()
            header_list = header.strip().split(separator)
            with open(output_file, "w") as out_fd:
                out_fd.write(header)
                for line in in_fd:
                    tmp_line = line.strip().split(separator)
                    data = np.array(map(float, tmp_line[1:]))
                    max_level = max(data)
                    if max_level == 0:
                        zero_max_lvl_list.append(tmp_line[0])

                        if verbose:
                            print("Zero max level for %s...Skipping..." % tmp_line[0])
                        continue

                    data /= max_level
                    output_string = tmp_line[0] + "\t"
                    output_string += "\t".join(map(str, data))
                    output_string += "\n"
                    out_fd.write(output_string)

        zero_max_lvl_list.write(zero_max_lvl_list_file)

    @staticmethod
    def extract_counts_by_max_level(input_file, output_prefix,
                                   separator="\t", verbose=True):
        output_file = "%s.divided_by_maxlvl" % output_prefix
        zero_max_lvl_list = IdList()

        zero_max_lvl_list_file = "%s.zero_max_lvl.ids" % output_prefix

        with open(input_file, "r") as in_fd:
            header = in_fd.readline()
            header_list = header.strip().split(separator)
            with open(output_file, "w") as out_fd:
                out_fd.write(header)
                for line in in_fd:
                    tmp_line = line.strip().split(separator)
                    data = np.array(map(float, tmp_line[1:]))
                    max_level = max(data)
                    if max_level == 0:
                        zero_max_lvl_list.append(tmp_line[0])

                        if verbose:
                            print("Zero max level for %s...Skipping..." % tmp_line[0])
                        continue

                    data /= max_level
                    output_string = tmp_line[0] + "\t"
                    output_string += "\t".join(map(str, data))
                    output_string += "\n"
                    out_fd.write(output_string)

        zero_max_lvl_list.write(zero_max_lvl_list_file)

    @staticmethod
    def divide_counts_by_several_base_level(input_file, output_prefix, base_levels,
                                            separator="\t", verbose=True,
                                            max_ratio_to_base_lvl=0.5):
        output_file = "%s.divided_by_max_baselvl" % output_prefix
        max_ratio_to_base_lvl_file = "%s.divided_by_max_baselvl.max_%f_ratio" % (output_prefix, max_ratio_to_base_lvl)
        zero_max_base_lvl_list = IdList()
        zero_max_base_lvl_list_file = "%s.zero_base_lvls.ids" % output_prefix
        max_ratio_to_base_lvl_fd = open(max_ratio_to_base_lvl_file, "w")
        with open(input_file, "r") as in_fd:
            header = in_fd.readline()
            header_list = header.strip().split(separator)

            data_base_lvl_index_list = []
            base_level_list = [base_levels] if isinstance(base_levels, str) else base_levels
            for level in base_level_list:
                data_base_lvl_index_list.append(header_list.index(level) - 1)

            with open(output_file, "w") as out_fd:
                out_fd.write(header)
                max_ratio_to_base_lvl_fd.write(header)
                for line in in_fd:
                    tmp_line = line.strip().split(separator)
                    data = np.array(map(float, tmp_line[1:]))
                    max_base_lvl = max(np.take(data, data_base_lvl_index_list))
                    if max_base_lvl == 0:
                        zero_max_base_lvl_list.append(tmp_line[0])
                        if verbose:
                            print("Zero max base level(s) for %s...Skipping..." % tmp_line[0])
                        continue

                    data /= max_base_lvl
                    output_string = tmp_line[0] + "\t"
                    output_string += "\t".join(map(str, data))
                    output_string += "\n"
                    if max(np.delete(data, data_base_lvl_index_list)) <= max_ratio_to_base_lvl:
                        max_ratio_to_base_lvl_fd.write(output_string)
                    out_fd.write(output_string)

        zero_max_base_lvl_list.write(zero_max_base_lvl_list_file)
        max_ratio_to_base_lvl_fd.close()
