#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import re
import argparse
from RouToolPa.Routines.File import split_filename, check_path



def read_entry(in_fd):

    name = in_fd.readline().strip()
    if name == "":
        return None, None, None, None
    sequence = in_fd.readline().strip()
    separator = in_fd.readline().strip()
    quality = in_fd.readline().strip()

    return name, sequence, separator, quality


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_se", action="store", dest="input_se",
                    help="Input file with SE reads")
parser.add_argument("-l", "--input_left", action="store", dest="input_left",
                    help="Input file with left reads")
parser.add_argument("-r", "--input_right", action="store", dest="input_right",
                    help="Input file with right reads")
parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", default="./",
                    help="Directory to write output")
parser.add_argument("-m", "--min_len", action="store", dest="min_len", type=int, default=1,
                    help="Minimum length of read to output")

args = parser.parse_args()
n_regexp = re.compile("N+$")

if args.input_se:

    se_directory, se_prefix, se_extension = split_filename(args.input_se)
    se_in_fd = open(args.input_se, "r")
    se_out_file = "%s%s.filtered%s" % (check_path(args.out_dir), se_prefix, se_extension)
    se_out_fd = open(se_out_file, "w")

    while True:
        name, sequence, separator, quality = read_entry(se_in_fd)
        if name is None:
            break
        match = n_regexp.search(sequence)
        if match is None:
            se_out_fd.write("%s\n%s\n%s\n%s\n" % (name, sequence, separator, quality))
        elif match.start() >= args.min_len:
            se_out_fd.write("%s\n%s\n%s\n%s\n" % (name, sequence[:match.start()+1], separator, quality[:match.start()+1]))
        else:
            continue

    se_in_fd.close()
    se_out_fd.close()

elif args.input_left and args.input_right:

    left_directory, left_prefix, left_extension = split_filename(args.input_left)
    right_directory, right_prefix, right_extension = split_filename(args.input_right)
    left_in_fd = open(args.input_left, "r")
    right_in_fd = open(args.input_right, "r")
    left_out_file = "%s%s.filtered_1%s" % (check_path(args.out_dir), left_prefix, left_extension)
    left_out_fd = open(left_out_file, "w")
    right_out_file = "%s%s.filtered_2%s" % (check_path(args.out_dir), right_prefix, right_extension)
    right_out_fd = open(right_out_file, "w")

    left_out_se_file = "%s%s.se.filtered_1%s" % (check_path(args.out_dir), left_prefix, left_extension)
    left_out_se_fd = open(left_out_se_file, "w")
    right_out_se_file = "%s%s.se.filtered_2%s" % (check_path(args.out_dir), right_prefix, right_extension)
    right_out_se_fd = open(right_out_se_file, "w")

    while True:
        left_name, left_sequence, left_separator, left_quality = read_entry(left_in_fd)
        right_name, right_sequence, right_separator, right_quality = read_entry(right_in_fd)
        if (left_name is None) or (right_name is None):
            break
        left_match = n_regexp.search(left_sequence)
        right_match = n_regexp.search(right_sequence)

        if (((left_match is None) and (len(left_sequence) >= args.min_len)) or ((left_match is not None) and (left_match.start() >= args.min_len))) and (((right_match is None) and (len(right_sequence) >= args.min_len)) or ((right_match is not None) and (right_match.start() >= args.min_len))):
            #print(left_match)
            #try:
            #    print(left_match.start())
            #except:
            #    pass
            #print(left_name)
            #print(left_sequence)
            #print(type(left_sequence if left_match is None else left_sequence[:left_match.start()]))
            #print(left_separator)
            #print(left_quality)
            #print(left_match.start())
            #print(left_quality if left_match is None else left_quality[:left_match.start()])
            left_out_fd.write("%s\n%s\n%s\n%s\n" % (left_name,
                                                     left_sequence if (left_match is None) else left_sequence[:left_match.start()],
                                                     left_separator,
                                                     left_quality if (left_match is None) else left_quality[:left_match.start()]))
            right_out_fd.write("%s\n%s\n%s\n%s\n" % (right_name,
                                                      right_sequence if right_match is None else right_sequence[:right_match.start()],
                                                      right_separator,
                                                      right_quality if right_match is None else right_quality[:right_match.start()]))
        elif ((left_match is None) and (len(left_sequence) >= args.min_len)) or ((left_match is not None) and ((left_match.start() + 1) >= args.min_len)):
            left_out_se_fd.write("%s\n%s\n%s\n%s\n" % (left_name,
                                                        left_sequence if left_match is None else left_sequence[:left_match.start()],
                                                        left_separator,
                                                        left_quality if left_match is None else left_quality[:left_match.start()]))
        elif ((right_match is None) and (len(right_sequence) >= args.min_len)) or ((right_match is not None) and ((right_match.start() + 1) >= args.min_len)):
            right_out_se_fd.write("%s\n%s\n%s\n%s\n" % (right_name,
                                                         right_sequence if right_match is None else right_sequence[:right_match.start()],
                                                         right_separator,
                                                         right_quality if right_match is None else right_quality[:right_match.start()]))
        else:
            continue
    left_in_fd.close()
    right_in_fd.close()
    left_out_fd.close()
    right_out_fd.close()
    left_out_se_fd.close()
    right_out_se_fd.close()

else:
    raise IOError("Wrong input")

