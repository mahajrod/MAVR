#!/usr/bin/env python
import os

from Tools.Abstract import Tool
from collections import OrderedDict


class Gblocks(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "Gblocks", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(input_type="codon",
                      min_seq_number_for_conserved_position=None,
                      min_seq_number_for_flank_position=None,
                      max_pos_number_for_noncons_contig_pos=None,
                      min_block_len=None,
                      allow_gaps="half",
                      save_postscript=False,
                      output_type="htm",
                      concatenate_blocks_from_aignments=None,
                      ):
        options = " -t=%s" % ("p" if input_type=='protein' else "d" if input_type=='dna' else "c")
        options += " -b1=%i" % min_seq_number_for_conserved_position if min_seq_number_for_conserved_position else ""
        options += " -b2=%i" % min_seq_number_for_flank_position if min_seq_number_for_flank_position else ""
        options += " -b3=%i" % max_pos_number_for_noncons_contig_pos if max_pos_number_for_noncons_contig_pos else ""
        options += " -b4=%i" % min_block_len if min_block_len else ""
        options += " -b5=%s" % ("h" if allow_gaps == 'half' else "a" if allow_gaps else "n")
        options += " -d=y" if save_postscript else ""
        options += " -a=y" if concatenate_blocks_from_aignments else ""
        options += " -a=y" if concatenate_blocks_from_aignments else ""
        options += " -p=%s" % ("y" if output_type == "htm" else "t" if output_type == "text" else "s" if output_type == "short_text" else "n")

        return options

    @staticmethod
    def convert_output_to_fasta(input_file, output_file):
        sed_string = "sed '/^[^>]/{s/\\ //g}' %s > %s" % (input_file, output_file)
        os.system(sed_string)

    @staticmethod
    def extract_block_coordinates(htm_file):
        with open(htm_file, "r") as in_fd:
            for line in in_fd:
                if line[:7] == "Flanks:":
                    tmp = line.split(":")[1].strip().split("  ")
                    block_coordinates = [] # 1-based
                    for i in range(0, len(tmp)/2):
                        #print tmp
                        block_coordinates.append((int(tmp[2*i][1:]), int(tmp[2*i+1][:-1])))

        return block_coordinates

    def parallel_run(self, input_dir, output_dir, output_prefix,
                     input_type="codon",
                     min_seq_number_for_conserved_position=None,
                     min_seq_number_for_flank_position=None,
                     max_pos_number_for_noncons_contig_pos=None,
                     min_block_len=None,
                     allow_gaps="half",
                     save_postscript=True,
                     output_type="htm",
                     threads=None,
                     ):

        if threads:
            self.threads = threads

        data_dir = "%s/data/" % output_dir
        postscript_dir = "%s/ps/" % output_dir
        results_dir = "%s/results/" % output_dir
        htm_dir = "%s/htm/" % output_dir

        for directory in output_dir, data_dir, postscript_dir, results_dir, htm_dir:
            self.safe_mkdir(directory)

        #input_files_list = map(os.path.abspath, self.make_list_of_path_to_files(input_directory))

        input_files_list = self.make_list_of_path_to_files(input_dir, return_absolute_paths=True)

        for entry in input_files_list:
            directory, prefix, extension = self.split_filename(entry)
            os.system("ln -s %s %s/%s.%s" % (entry, data_dir, prefix, extension))

        data_files_list = self.make_list_of_path_to_files(input_dir, return_absolute_paths=True)

        common_options = self.parse_options(input_type=input_type,
                                            min_seq_number_for_conserved_position=min_seq_number_for_conserved_position,
                                            min_seq_number_for_flank_position=min_seq_number_for_flank_position,
                                            max_pos_number_for_noncons_contig_pos=max_pos_number_for_noncons_contig_pos,
                                            min_block_len=min_block_len,
                                            allow_gaps=allow_gaps,
                                            save_postscript=save_postscript,
                                            output_type=output_type,
                                            concatenate_blocks_from_aignments=None)
        options_list = []
        for data_file in data_files_list:
            options = " %s" % data_file
            options += " %s" % common_options
            options_list.append(options)

        self.parallel_execute(options_list=options_list)

        block_coordinates = OrderedDict()

        for filename in data_files_list:
            data_dir, prefix, extension = self.split_filename(filename)
            blocks_file = "%s-gb" % filename
            htm_file = "%s-gb.htm" % filename
            postscript_file = "%s-gbPS" % filename

            block_coordinates[prefix] = self.extract_block_coordinates(htm_file)
            os.system("mv %s %s/%s.ps" % (postscript_file, postscript_dir, prefix))
            os.system("mv %s %s/%s.htm" % (htm_file, htm_dir, prefix))
            self.convert_output_to_fasta(blocks_file, "%s/%s%s" % (results_dir, prefix, extension))

        block_coordinates_file = "%s.block.coordinates" % output_prefix

        with open(block_coordinates_file, "w") as block_fd:
            for entry in block_coordinates:
                coordinates_string = ";".join(map(lambda s: "%i,%i" % (s[0], s[1]), block_coordinates[entry]))
                block_fd.write("%s\t%s\n" % (entry, coordinates_string))
