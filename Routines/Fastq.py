__author__ = 'mahajrod'
import os
import re
import sys
import pickle

from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from CustomCollections.GeneralCollections import TwoLvlDict, SynDict, IdList, IdSet
from Routines.Functions import output_dict
from Routines.File import FileRoutines


class FastQRoutines(FileRoutines):

    def __init__(self):
        pass

    @staticmethod
    def reverse_complement(in_file, out_file):
        with open(in_file, "r") as in_fd:
            with open(out_file, "w") as out_fd:
                for line in in_fd:
                    out_fd.write(line)
                    out_fd.write(str(Seq(in_fd.next().strip()).reverse_complement()))
                    out_fd.write("\n")
                    out_fd.write(in_fd.next())
                    out_fd.write(in_fd.next().strip()[::-1])
                    out_fd.write("\n")

    @staticmethod
    def parse_illumina_name(string):
        name_list = string[1:].split(" ")

    def remove_tiles_from_fastq(self, forward_reads, black_list_forward_tiles_list,
                                reverse_reads, black_list_reverse_tiles_list, output_prefix):

        filtered_paired_forward_pe = "%s.ok.pe_1.fastq" % output_prefix
        filtered_forward_se = "%s.ok.forward.se.fastq" % output_prefix
        filtered_out_forward_se = "%s.bad.forward.fastq" % output_prefix

        filtered_paired_reverse_pe = "%s.ok.pe_2.fastq" % output_prefix
        filtered_reverse_se = "%s.ok.reverse.se.fastq" % output_prefix
        filtered_out_reverse_se = "%s.bad.reverse.fastq" % output_prefix

        forward_input_fd = self.metaopen(forward_reads, "r")
        reverse_input_fd = self.metaopen(reverse_reads, "r")

        filtered_paired_forward_pe_fd = self.metaopen(filtered_paired_forward_pe, "w")
        filtered_forward_se_fd = self.metaopen(filtered_forward_se, "w")
        filtered_out_forward_se_fd = self.metaopen(filtered_out_forward_se, "w")

        filtered_paired_reverse_pe_fd = self.metaopen(filtered_paired_reverse_pe, "w")
        filtered_reverse_se_fd = self.metaopen(filtered_reverse_se, "w")
        filtered_out_reverse_se_fd = self.metaopen(filtered_out_reverse_se, "w")

        for line in forward_input_fd:
            name_list = line.strip()[1:].split(":")
            #instrument_id = name_list[0]
            #run_id = name_list[1]
            #flowcell_id = name_list[2]
            #lane_id = name_list[3]
            tile_id = name_list[4]

            if (tile_id in black_list_forward_tiles_list) and (tile_id in black_list_reverse_tiles_list):
                filtered_out_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_out_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_out_reverse_se_fd.write(reverse_input_fd.next())

            elif (tile_id in black_list_forward_tiles_list) and (not(tile_id in black_list_reverse_tiles_list)):
                filtered_out_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_out_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_reverse_se_fd.write(reverse_input_fd.next())

            elif (not (tile_id in black_list_forward_tiles_list)) and (tile_id in black_list_reverse_tiles_list):
                filtered_forward_se_fd.write(line)
                for i in range(0, 3):
                    filtered_forward_se_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_out_reverse_se_fd.write(reverse_input_fd.next())

            else:
                filtered_paired_forward_pe_fd.write(line)
                for i in range(0, 3):
                    filtered_paired_forward_pe_fd.write(forward_input_fd.next())
                for i in range(0, 4):
                    filtered_paired_reverse_pe_fd.write(reverse_input_fd.next())

        filtered_paired_forward_pe_fd.close()
        filtered_forward_se_fd.close()
        filtered_out_forward_se_fd.close()

        filtered_paired_reverse_pe_fd.close()
        filtered_reverse_se_fd.close()
        filtered_out_reverse_se_fd.close()

    def combine_fastq_files(self, samples_directory, sample, output_directory,
                            use_links_if_merge_not_necessary=True):
        sample_dir = "%s/%s/" % (samples_directory, sample)

        filetypes, forward_files, reverse_files, se_files = self.make_lists_forward_and_reverse_files(sample_dir)

        uncompresed = True
        if len(filetypes) == 1:
            if ("fq.gz" in filetypes) or ("fastq.gz" in filetypes):
                command = "zcat"
                uncompresed = False
            elif ("fq.bz2" in filetypes) or ("fastq.bz2" in filetypes):
                command = "bzcat"
                uncompresed = False
            else:
                command = "cat"

            merged_forward = "%s/%s_1.fq" % (output_directory, sample)
            merged_reverse = "%s/%s_2.fq" % (output_directory, sample)
            merged_se = "%s/%s.se.fq" % (output_directory, sample)

            if use_links_if_merge_not_necessary and (len(forward_files) == 1) and (len(reverse_files) == 1) and (uncompresed == True):
                os.system("ln -s %s %s" % (forward_files[0], merged_forward))
                os.system("ln -s %s %s" % (reverse_files[0], merged_reverse))
                if len(se_files) == 1:
                    os.system("ln -s %s %s" % (se_files[0], merged_se))
                    return merged_forward, merged_reverse, merged_se
                elif len(se_files) > 0:
                    os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                    return merged_forward, merged_reverse, None
            else:
                if (len(forward_files) > 0) and (len(reverse_files) > 0):
                    os.system("%s %s > %s" % (command, " ".join(forward_files), merged_forward))
                    os.system("%s %s > %s" % (command, " ".join(reverse_files), merged_reverse))
                    if len(se_files) > 0:
                        os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                        return merged_forward, merged_reverse, merged_se
                    else:
                        return merged_forward, merged_reverse, None
                if len(se_files) > 0:
                    os.system("%s %s > %s" % (command, " ".join(se_files), merged_se))
                    return None, None, merged_se
                else:
                    raise IOError("No input files were found")
        else:
            raise IOError("Extracting from mix of archives in not implemented yet")

    @staticmethod
    def filter_se_by_length(input_file, filtered_file, filtered_out_file, min_len=None, max_len=None):

        if min_len and max_len:
            def expression(line):
                return min_len <= len(line) <= max_len
        elif min_len:
            def expression(line):
                return min_len <= len(line)
        elif max_len:
            def expression(line):
                return len(line) <= max_len
        else:
            raise ValueError("Both minimum and maximum thresholds for read length were not set")

        with open(input_file, "r") as in_fd:
            with open(filtered_file, "w") as filtered_fd:
                with open(filtered_out_file, "w") as filtered_out_fd:
                    for read_name in in_fd:
                        read = in_fd.next()
                        #print len(read.strip())
                        #print expression(read.strip())
                        if expression(read.strip()):
                            filtered_fd.write(read_name)
                            filtered_fd.write(read)
                            filtered_fd.write(in_fd.next())
                            filtered_fd.write(in_fd.next())
                        else:
                            filtered_out_fd.write(read_name)
                            filtered_out_fd.write(read)
                            filtered_out_fd.write(in_fd.next())
                            filtered_out_fd.write(in_fd.next())
