#!/usr/bin/env python
__author__ = 'mahajrod'


def read_ids(filename, header=False):
    #reads ids from file with one id per line
    id_list = []
    with open(filename, "r") as in_fd:
        if header:
            header_str = in_fd.readline().strip()
        for line in in_fd:
            id_list.append(line.strip())
    return id_list


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
