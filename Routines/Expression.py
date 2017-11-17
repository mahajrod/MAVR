#!/usr/bin/env python
__author__ = 'mahajrod'
import numpy as np
from CustomCollections.GeneralCollections import IdList, SynDict


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

    @staticmethod
    def calculate_fpkm_for_count_table(count_table_file, transcript_length_file, output_file,
                                       separator="\t"):
        length_dict = SynDict(filename=transcript_length_file, expression=int, comments_prefix="#")

        with open(count_table_file, "r") as in_fd:
            header_list = in_fd.readline().strip().split(separator)

            samples_list = header_list[1:]
            gene_list = IdList()
            count_list = []
            for line in in_fd:
                tmp = line.strip().split(separator)
                gene_list.append(tmp[0])
                count_list.append(map(float, tmp[1:]))

            per_sample_total_counts = []

            for sample_index in range(0, len(samples_list)):
                total_counts = 0
                for gene_index in range(0, len(count_list)):
                    total_counts += count_list[gene_index][sample_index]
                per_sample_total_counts.append(total_counts)

        with open(output_file, "w") as out_fd:
            out_fd.write(separator.join(header_list) + "\n")
            for gene_index in range(0, len(count_list)):
                normalized_counts_list = []
                for sample_index in range(0, len(samples_list)):
                    gene_count = count_list[gene_index][sample_index] * (10**9) / length_dict[gene_list[gene_index]] / per_sample_total_counts[sample_index]
                    normalized_counts_list.append(gene_count)
                out_fd.write("%s\t%s\n" % (gene_list[gene_index], "\t".join(map(str, normalized_counts_list))))
