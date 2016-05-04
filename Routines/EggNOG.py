#!/usr/bin/env python
import os

from Bio import SeqIO

from Routines import AlignmentRoutines, FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict


class EggNOGRoutines:
    def __init__(self):

        pass

    @staticmethod
    def edit_profile_names_in_fam_file(input_fam_file, output_fam_file, comments_prefix="#"):
        if comments_prefix:
            comments_prefix_len = len(comments_prefix)

        with open(input_fam_file, "r") as in_fd:
            with open(output_fam_file, "w") as out_fd:
                for line in in_fd:
                    if comments_prefix:
                        if line[:comments_prefix_len] == comments_prefix:
                            out_fd.write(line)
                    line_list = line.split("\t")
                    fam_name = line_list[0].split(".")[1]
                    out_fd.write("%s\t%s" % (fam_name, line_list[1]))

    @staticmethod
    def convert_members_tsv_to_fam(input_file, output_file):
        cmd = "awk -F'\t' '{printf \"%%s\\t%%s\\n\",$2,$6 }' %s > %s" % (input_file, output_file)
        os.system(cmd)

    @staticmethod
    def extract_proteins_from_alignments(dir_with_alignments, output_dir):
        input_files = FileRoutines.make_list_of_path_to_files([dir_with_alignments] if isinstance(dir_with_alignments, str) else dir_with_alignments)
        out_dir = FileRoutines.check_path(output_dir)
        FileRoutines.save_mkdir(output_dir)
        for filename in input_files:
            filename_list = FileRoutines.split_filename(filename)
            output_file = "%s%s%s" % (out_dir, filename_list[1], filename_list[2])
            AlignmentRoutines.extract_sequences_from_alignment(filename, output_file)



