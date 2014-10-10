__author__ = 'mahajrod'

import sys


def check_path(path_to_check):
    #returns path with / at end or blank path
    if path_to_check != "":
        if path_to_check[-1] != "/":
            return path_to_check + "/"
    return path_to_check


def output_dict(dictionary, out_file="out.t", write=False):

    out_fd = open(out_file, "w") if write else sys.stdout

    for key in dictionary:
        out_fd.write("%s\t%s\n" % (key, str(dictionary[key])))

    if out_fd is not sys.stdout:
        out_fd.close()