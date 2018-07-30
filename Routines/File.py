#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import sys
import bz2
import gzip
from collections import Iterable, OrderedDict

from CustomCollections.GeneralCollections import IdSet,  IdList, SynDict


class FileRoutines:
    def __init__(self):
        self.filetypes_dict = {"fasta": [".fa", ".fasta", ".fa", ".pep", ".cds"],
                               "fastq": [".fastq", ".fq"],
                               "genbank": [".gb", ".genbank"],
                               "newick": [".nwk"],
                               "gz": [".gz"],
                               "bzip": [".bz2"]}

    @staticmethod
    def metaopen(filename, flags):
        if filename[-3:] == ".gz":
            return gzip.open(filename, flags)
        elif filename[-4:] == ".bz2":
            return bz2.open(filename, flags)
        else:
            return open(filename, flags)

    @staticmethod
    def safe_mkdir(dirname, 
                   description_filename=None, description_text=None,
                   readme_filename=None, readme_text=None):
        try:
            os.mkdir(dirname)
        except OSError:
            pass
        
        if not(description_filename is None):
            description_filename = "%s/%s" % (dirname, description_filename)
    
            if not os.path.isfile(description_filename):
                with open(description_filename, "w") as descr_fd:
                    if not (description_text is None):
                        descr_fd.write(description_text)
                        
        if not(readme_filename is None):
            readme_filename = "%s/%s" % (dirname, readme_filename)
    
            if not os.path.isfile(readme_filename):
                with open(readme_filename, "w") as descr_fd:
                    if not (readme_text is None):
                        descr_fd.write(readme_text)

    def recursive_mkdir(self, dir_dict, out_dir=None,
                        description_filename=None, description_text=None,
                        readme_filename=None, readme_text=None):
        if not(out_dir is None):
            self.safe_mkdir(out_dir)
        for directory in dir_dict:
            dirname = directory if out_dir is None else "%s/%s" % (out_dir, directory)
            self.safe_mkdir(dirname,
                            description_filename=description_filename,
                            description_text=description_text,
                            readme_filename=readme_filename,
                            readme_text=readme_text)

            if isinstance(dir_dict[directory], dict):
                self.recursive_mkdir(dir_dict[directory],
                                     out_dir=dirname,
                                     description_filename=description_filename,
                                     description_text=description_text,
                                     readme_filename=readme_filename,
                                     readme_text=readme_text)



    def detect_filetype_by_extension(self, filename, filetypes_dict=None):
        filetypes = filetypes_dict if filetypes_dict else self.filetypes_dict
        directory, prefix, extension = split_filename(filename)
        for filetype in filetypes:
            if extension in filetypes[filetype]:
                return filetype
        return None

    @staticmethod
    def check_path(path_to_check):
        #print (path_to_check)
        #returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    @staticmethod
    def check_dir_path(path_to_check):
        #print (path_to_check)
        #returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    @staticmethod
    def split_filename(filepath):
        directory, basename = os.path.split(filepath)
        prefix, extension = os.path.splitext(basename)
        return directory, prefix, extension

    def make_list_of_path_to_files(self, list_of_dirs_and_files, expression=None, recursive=False,
                                   return_absolute_paths=True):
        file_list = []
        for entry in [list_of_dirs_and_files] if isinstance(list_of_dirs_and_files, str) else list_of_dirs_and_files:
            if os.path.isdir(entry):
                files_in_dir = ["%s%s" % (self.check_path(entry), filename)
                                for filename in sorted(filter(expression, os.listdir(entry))
                                                       if expression else os.listdir(entry))]
                if recursive:
                    for filename in files_in_dir:
                        if os.path.isdir(filename):
                            file_list += self.make_list_of_path_to_files([filename],
                                                                         expression=expression,
                                                                         recursive=recursive)
                        else:
                            file_list.append(filename)
                else:
                    file_list += files_in_dir
            elif os.path.exists(entry):
                if expression:
                    if expression(os.path.abspath(entry)):
                        file_list.append(os.path.abspath(entry))
                        print "BBBBBBBBBBBBBB"
                        print entry
                    else:
                        print "AAAAAAAAAAA"
                        print entry
                else:
                    file_list.append(os.path.abspath(entry))
            else:
                print("%s does not exist" % entry)

        return map(os.path.abspath, file_list) if return_absolute_paths else file_list

    def make_list_of_path_to_files_from_string(self, input_string, file_separator=",",
                                               expression=None, recursive=False):
        return self.make_list_of_path_to_files(input_string.split(file_separator),
                                               expression=expression,
                                               recursive=recursive)

    @staticmethod
    def check_extension(filename, extension_list=[]):
        if extension_list:
            for extension in extension_list:
                if extension == filename[-len(extension):]:
                    return True
        else:
            return True

        return False

    def make_list_of_path_to_files_by_extension(self, list_of_dirs_and_files, extension_list=[], recursive=False,
                                                return_absolute_paths=True):

        def check_extension(filename):
            return self.check_extension(filename, extension_list)

        return self.make_list_of_path_to_files(list_of_dirs_and_files, expression=check_extension,
                                               recursive=recursive, return_absolute_paths=return_absolute_paths)

    @staticmethod
    def read_synonyms_dict(filename, header=False, separator="\t",
                           split_values=False, values_separator=",", key_index=0, value_index=1):
        # reads synonyms from file
        synonyms_dict = OrderedDict()
        with open(filename, "r") as in_fd:
            if header:
                header_str = in_fd.readline().strip()
            for line in in_fd:
                tmp = line.strip().split(separator)
                #print line
                key, value = tmp[key_index], tmp[value_index]
                if split_values:
                    value = value.split(values_separator)
                synonyms_dict[key] = value
        return synonyms_dict

    @staticmethod
    def read_ids(filename, header=False, close_after_if_file_object=False):
        #reads ids from file or file object with one id per line
        id_list = []

        in_fd = filename if isinstance(filename, file) else open(filename, "r")

        if header:
            header_str = in_fd.readline().strip()
        for line in in_fd:
            id_list.append(line.strip())
        if (not isinstance(filename, file)) or close_after_if_file_object:
            in_fd.close()

        return id_list

    @staticmethod
    def read_tsv_as_rows_list(filename, header=False, header_prefix="#", separator="\t"):
        tsv_list = []
        with open(filename, "r") as tsv_fd:
            if header:
                header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
            for line in tsv_fd:
                tsv_list.append(line.strip().split(separator))
        return header_list, tsv_list if header else tsv_list

    @staticmethod
    def read_tsv_as_columns_dict(filename, header_prefix="#", separator="\t"):
        tsv_dict = OrderedDict()
        with open(filename, "r") as tsv_fd:
            header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
            number_of_columns = len(header_list)
            for column_name in header_list:
                tsv_dict[column_name] = []
            for line in tsv_fd:
                tmp_line = line.strip().split(separator)
                for i in range(0, number_of_columns):
                    tsv_dict[header_list[i]].append(tmp_line[i])
        return tsv_dict

    @staticmethod
    def tsv_split_by_column(tsv_file, column_number, separator="\t", header=False, outfile_prefix=None):
        # column number should start from 0
        header_string = None
        splited_name = tsv_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd_dict = {}

        with open(tsv_file, "r") as in_fd:
            if header:
                header_string = in_fd.readline()
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in out_fd_dict:
                    print (line_str[column_number])
                    suffix = line_str[column_number].replace(" ", "_")
                    out_name = "%s_%s.%s" % (out_prefix, suffix, extension)
                    out_fd_dict[line_str[column_number]] = open(out_name, "w")
                    if header:
                        out_fd_dict[line_str[column_number]].write(header_string)
                out_fd_dict[line_str[column_number]].write(line)

        for entry in out_fd_dict:
            out_fd_dict[entry].close()

    @staticmethod
    def tsv_extract_by_column_value(tsv_file, column_number, column_value, separator="\t",
                                    header=False, outfile_prefix=None):
        # column number should start from 0
        # column_value should be string or list of strings

        splited_name = tsv_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd = open("%s.%s" % (out_prefix, extension), "w")

        column_value_list = column_value if isinstance(column_value, Iterable) else []
        with open(tsv_file, "r") as in_fd:
            if header:
                out_fd.write(in_fd.readline())
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] in column_value_list:
                    out_fd.write(line)
        out_fd.close()

    @staticmethod
    def tsv_remove_by_column_value(tsv_file, column_number, column_value, separator="\t",
                                   header=False, outfile_prefix=None):
        # column number should start from 0
        # column_value should be string or list of strings

        splited_name = tsv_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd = open("%s.%s" % (out_prefix, extension), "w")

        column_value_list = column_value if isinstance(column_value, Iterable) else []
        with open(tsv_file, "r") as in_fd:
            if header:
                out_fd.write(in_fd.readline())
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in column_value_list:
                    out_fd.write(line)
        out_fd.close()

    @staticmethod
    def make_lists_forward_and_reverse_files(sample_dir, filename_fragment_to_mark_se_reads=".se.", input_is_se=False):
        file_list = sorted(os.listdir(sample_dir))
        filtered_filelist = []
        filetypes = set()
        for entry in file_list:
            if entry[-3:] == ".fq" or entry[-6:] == ".fastq":
                filetypes.add("fq")
            elif entry[-6:] == ".fq.gz" or entry[-9:] == ".fastq.gz":
                filetypes.add("fq.gz")
            elif entry[-7:] == ".fq.bz2" or entry[-10:] == ".fastq.bz2":
                filetypes.add("fq.bz2")
            else:
                continue
            filtered_filelist.append("%s/%s" % (sample_dir, entry))

        if input_is_se:
            return filetypes, [], [], filtered_filelist

        single_end_filelist = []
        paired_end_filelist = []

        for entry in filtered_filelist:
            if filename_fragment_to_mark_se_reads in entry:
                single_end_filelist.append(entry)
            else:
                paired_end_filelist.append(entry)

        forward_files = paired_end_filelist[::2]
        reverse_files = paired_end_filelist[1:][::2]

        if len(filetypes) > 1:
            print("WARNING: mix of archives of different types and/or uncompressed files")

        return filetypes, forward_files, reverse_files, single_end_filelist

    @staticmethod
    def get_sample_list(samples_directory):
        samples = sorted(os.listdir(samples_directory))
        sample_list = []
        for sample in samples:
            if os.path.isdir("%s/%s" % (samples_directory, sample)):
                sample_list.append(sample)
        return sample_list

    @staticmethod
    def extract_ids_from_file(input_file, output_file=None, header=False, column_separator="\t",
                              comments_prefix="#", column_number=None):
        id_list = IdList()
        id_list.read(input_file, column_separator=column_separator, comments_prefix=comments_prefix,
                     column_number=column_number, header=header)
        if output_file:
            id_list.write(output_file, header=header)
        return id_list

    @staticmethod
    def intersect_ids_from_files(files_with_ids_from_group_a, files_with_ids_from_group_b,
                                 result_file=None, mode="common"):
        a = IdSet()
        b = IdSet()

        if mode == "common":
            expression = lambda a, b: a & b
        elif mode == "only_a":
            expression = lambda a, b: a - b
        elif mode == "only_b":
            expression = lambda a, b: b - a
        elif mode == "not_common":
            expression = lambda a, b: a ^ b
        elif mode == "combine":
            expression = lambda a, b: a | b

        #print(files_with_ids_from_group_a)
        for filename in [files_with_ids_from_group_a] if isinstance(files_with_ids_from_group_a, str) else files_with_ids_from_group_a:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#")
            a = a | id_set

        for filename in [files_with_ids_from_group_b] if isinstance(files_with_ids_from_group_b, str) else files_with_ids_from_group_b:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#")
            b = b | id_set

        result_fd = open(result_file, "w") if result_file else sys.stdout
        if mode != "count":
            final_set = IdSet(expression(a, b))
            final_set.write(result_fd)
        else:
            result_fd.write("Group_A\t%i\nGroup_B\t%i\nCommon\t%i\nOnly_group_A\t%i\nOnly_group_B\t%i\nNot_common\t%i\nAll\t%i\n" %
                            (len(a), len(b), len(a & b), len(a - b), len(b - a), len(a ^ b), len(a | b)))

    @staticmethod
    def intersect_ids(list_of_group_a, list_of_group_b, mode="common"):
        # possible modes: common, only_a, only_b, not_common,  combine, count
        a = IdSet()
        b = IdSet()

        if mode == "common":
            expression = lambda a, b: a & b
        elif mode == "only_a":
            expression = lambda a, b: a - b
        elif mode == "only_b":
            expression = lambda a, b: b - a
        elif mode == "not_common":
            expression = lambda a, b: a ^ b
        elif mode == "combine":
            expression = lambda a, b: a | b

        for id_list in list_of_group_a:
            a = a | IdSet(id_list)

        for id_list in list_of_group_b:
            b = b | IdSet(id_list)

        if mode != "count":
            return IdSet(expression(a, b))
        else:
            return len(a), len(b), len(a & b), len(a - b), len(b - a), len(a ^ b), len(a | b)

    @staticmethod
    def split_by_column(input_file, column_number, separator="\t", header=False, outfile_prefix=None,
                        use_column_value_as_prefix=False, sorted_input=False):
        # column number should start from 0
        # use sorted input to reduce number of simalteniously open files
        header_string = None
        splited_name = input_file.split(".")
        extension = splited_name[-1] if len(splited_name) > 1 else ""
        out_prefix = outfile_prefix if outfile_prefix is not None \
            else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
        out_fd_dict = {}
        previous_value = None
        with open(input_file, "r") as in_fd:
            if header:
                header_string = in_fd.readline()
            for line in in_fd:
                line_str = line.strip().split(separator)
                if line_str[column_number] not in out_fd_dict:
                    print (line_str[column_number])
                    if previous_value and sorted_input:
                        out_fd_dict[previous_value].close()
                    suffix = line_str[column_number].replace(" ", "_")
                    out_name = "%s.%s" % (suffix, extension) if use_column_value_as_prefix else "%s.%s.%s" % (out_prefix, suffix, extension)
                    out_fd_dict[line_str[column_number]] = open(out_name, "w")
                    if header:
                        out_fd_dict[line_str[column_number]].write(header_string)

                out_fd_dict[line_str[column_number]].write(line)
                previous_value = line_str[column_number]
        if sorted_input:
            out_fd_dict[previous_value].close()
        else:
            for entry in out_fd_dict:
                out_fd_dict[entry].close()

    @staticmethod
    def get_sample_list(sample_dir, sample_list=None):
        samples = []
        if sample_list:
            return [sample_list] if isinstance(sample_list, str) else sample_list
        else:
            dir_list = os.listdir(sample_dir)
            for directory in dir_list:
                if os.path.isdir("%s/%s" % (sample_dir, directory)):
                    samples.append(directory)

            return samples

    @staticmethod
    def replace_column_value_by_syn(input_file, syn_file, out_file, column=0, comment_prefix=None, separator="\t",
                                    syn_header=False, syn_separator="\t",
                                    syn_key_index=0, syn_value_index=1, syn_comment_prefix=None):
        syn_dict = SynDict(filename=syn_file, header=syn_header, separator=syn_separator, key_index=syn_key_index,
                           value_index=syn_value_index, comments_prefix=syn_comment_prefix)
        if comment_prefix:
            comment_prefix_len = len(comment_prefix)
        line_number = 0
        replaced = 0
        not_replaced = 0
        with open(input_file, "r") as in_fd:
            with open(out_file, "w") as out_fd:
                for line in in_fd:
                    line_number += 1
                    if comment_prefix:
                        if line[0:comment_prefix_len] == comment_prefix:
                            out_fd.write(line)
                            continue
                    line_list = line.strip("\n").split(separator)
                    if len(line_list) < column + 1:
                        sys.stderr.write("WARNING!!! Line %i doesn't have column %i\n" % (line_number, column))
                    if line_list[column] in syn_dict:
                        replaced += 1
                        line_list[column] = syn_dict[line_list[column]]
                    else:
                        not_replaced += 1

                    out_fd.write(separator.join(line_list))
                    out_fd.write("\n")

        sys.stderr.write("Replaced: %i\nNot replaced: %i\n" % (replaced, not_replaced))

    @staticmethod
    def combine_syn_dicts(list_of_syn_dict):
        combined_dict = SynDict()
        for syn_dict in list_of_syn_dict:
            for key in syn_dict:
                if key in combined_dict:
                    combined_dict[key] += syn_dict[key]
                else:
                    combined_dict[key] = syn_dict[key]

        return combined_dict

    def combine_syn_dicts_from_file(self, list_of_syndict_files, output, key_index=0, value_index=1, separator="\t",
                                    values_separator=",", header=False):
        list_of_syn_dicts = []
        for filename in list_of_syndict_files:
            list_of_syn_dicts.append(SynDict(filename=filename, key_index=key_index, value_index=value_index,
                                             separator=separator, values_separator=values_separator, split_values=True,
                                             allow_repeats_of_key=True, header=header))

        merged_dict = self.combine_syn_dicts(list_of_syn_dicts)

        merged_dict.write(output, splited_values=True)

    def add_add_new_column_by_key_column(self, table_file, syn_dict_file, key_column,  output_file, new_column_name=None,
                                         separator='\t', absent_value="."):
        column_syn_dict = SynDict(filename=syn_dict_file, allow_repeats_of_key=True, values_separator="@")
        with open(table_file, "r") as in_fd, open(output_file, "w") as out_fd:
            if new_column_name:
                header_line = in_fd.readline().strip() + "\t%s\n" % new_column_name
                out_fd.write(header_line)
                for line in in_fd:
                    line_list = line.strip().split(separator)
                    if line_list[key_column] in column_syn_dict:
                        print (line_list[key_column])
                        print (column_syn_dict[line_list[key_column]])
                    line_list.append(absent_value if line_list[key_column] not in column_syn_dict else "|".join(column_syn_dict[line_list[key_column]]))
                    out_fd.write(separator.join(line_list) + "\n")

    @staticmethod
    def label_column_in_file(input_file, label, column_index, output_file, column_separator="\t",
                             label_position="first", label_separator="@",
                             comments_prefix="#"):
        with open(input_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:
                for line in in_fd:
                    if line[0] == comments_prefix:
                        out_fd.write(line)
                        continue

                    line_list = line.strip("\n").split(column_separator)
                    if label_position == "first":
                        line_list[column_index] = "%s%s%s" % (label, label_separator, line_list[column_index])
                    elif label_position == "last":
                        line_list[column_index] = "%s%s%s" % (line_list[column_index], label_separator, label)
                    else:
                        raise ValueError("ERROR!!!Unrecognized label position")

                    out_fd.write(column_separator.join(line_list) + "\n")


filetypes_dict = {"fasta": [".fa", ".fasta", ".fa", ".pep", ".cds"],
                  "fastq": [".fastq", ".fq"],
                  "genbank": [".gb", ".genbank"],
                  "newick": [".nwk"]}


def safe_mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError:
        pass


def detect_filetype_by_extension(filename, filetypes_dict=filetypes_dict):
    directory, prefix, extension = split_filename(filename)
    for filetype in filetypes_dict:
        if extension in filetypes_dict[filetype]:
            return filetype
    return None


def check_path(path_to_check):
    #print (path_to_check)
    #returns path with / at end or blank path
    if path_to_check != "":
        if path_to_check[-1] != "/":
            return path_to_check + "/"
    return path_to_check


def split_filename(filepath):
    directory, basename = os.path.split(filepath)
    prefix, extension = os.path.splitext(basename)
    return directory, prefix, extension


def make_list_of_path_to_files(list_of_dirs_and_files, expression=None):

    paths_list = []
    for entry in list_of_dirs_and_files:
        #print entry
        if os.path.isdir(entry):
            files_in_dir = sorted(filter(expression, os.listdir(entry)) if expression else os.listdir(entry))
            for filename in files_in_dir:
                paths_list.append("%s%s" % (check_path(entry), filename))
        elif os.path.exists(entry):
            if expression:
                if expression(entry):
                    paths_list.append(entry)
            else:
                paths_list.append(entry)
        else:
            print("%s does not exist" % entry)

    return paths_list


def read_synonyms_dict(filename, header=False, separator="\t",
                       split_values=False, values_separator=",", key_index=0, value_index=1):
    # reads synonyms from file
    synonyms_dict = OrderedDict()
    with open(filename, "r") as in_fd:
        if header:
            header_str = in_fd.readline().strip()
        for line in in_fd:
            tmp = line.strip().split(separator)
            #print line
            key, value = tmp[key_index], tmp[value_index]
            if split_values:
                value = value.split(values_separator)
            synonyms_dict[key] = value
    return synonyms_dict


def read_ids(filename, header=False, close_after_if_file_object=False):
    #reads ids from file or file object with one id per line
    id_list = []

    in_fd = filename if isinstance(filename, file) else open(filename, "r")

    if header:
        header_str = in_fd.readline().strip()
    for line in in_fd:
        id_list.append(line.strip())
    if (not isinstance(filename, file)) or close_after_if_file_object:
        in_fd.close()

    return id_list


def read_tsv_as_rows_list(filename, header=False, header_prefix="#", separator="\t"):
    tsv_list = []
    with open(filename, "r") as tsv_fd:
        if header:
            header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
        for line in tsv_fd:
            tsv_list.append(line.strip().split(separator))
    return header_list, tsv_list if header else tsv_list


def read_tsv_as_columns_dict(filename, header_prefix="#", separator="\t"):
    tsv_dict = OrderedDict()
    with open(filename, "r") as tsv_fd:
        header_list = tsv_fd.readline().strip()[len(header_prefix):].split(separator)
        number_of_columns = len(header_list)
        for column_name in header_list:
            tsv_dict[column_name] = []
        for line in tsv_fd:
            tmp_line = line.strip().split(separator)
            for i in range(0, number_of_columns):
                tsv_dict[header_list[i]].append(tmp_line[i])
    return tsv_dict

def tsv_extract_by_column_value(tsv_file, column_number, column_value, separator="\t",
                                header=False, outfile_prefix=None):
    # column number should start from 0
    # column_value should be string or list of strings

    splited_name = tsv_file.split(".")
    extension = splited_name[-1] if len(splited_name) > 1 else ""
    out_prefix = outfile_prefix if outfile_prefix is not None \
        else ".".join(splited_name[:-1]) if len(splited_name) > 1 else splited_name[0]
    out_fd = open("%s.%s" % (out_prefix, extension), "w")

    column_value_list = column_value if isinstance(column_value, Iterable) else []
    with open(tsv_file, "r") as in_fd:
        if header:
            out_fd.write(in_fd.readline())
        for line in in_fd:
            line_str = line.strip().split(separator)
            if line_str[column_number] in column_value_list:
                out_fd.write(line)
    out_fd.close()
