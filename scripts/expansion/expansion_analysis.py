#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with blast alignment. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="gene_clusters.t",
                    help="Output file with genes clustered in families. Default: gene_clusters.t")
parser.add_argument("-a", "--solar_script_path", action="store", dest="solar_script_path", default="solar.pl",
                    help="Path to solar script. Default: solar.pl")
parser.add_argument("-g", "--hcluster_sg_path", action="store", dest="hcluster_sg_path", default="hcluster_sg",
                    help="Path to hcluster_sg binary. Default: hcluster_sg")
parser.add_argument("-l", "--solar_output_file", action="store", dest="solar_output_file", default="solar_output.t",
                    help="File with solar output. Default: solar_output.t")
parser.add_argument("-m", "--modified_solar_output_file", action="store", dest="modified_solar_output_file", default="modified_solar_output.t",
                    help="File with modified solar output. Default: modified_solar_output.t")
parser.add_argument("-b", "--bit_score_file", action="store", dest="bit_score_file", default="bit_score.t",
                    help="File with solar bitscores. Default: bit_score.t")
parser.add_argument("-e", "--self_alignment_score_file", action="store", dest="self_alignment_score_file", default="self_alignment_score.t",
                    help="File with self alignment scores. Default: self_alignment_score.t")
parser.add_argument("-t", "--temp_file", action="store", dest="temp_file", default="temp_file.t",
                    help="Temp file. Default: temp_file.t")
parser.add_argument("-c", "--hscore_file", action="store", dest="hscore_file", default="hscore_file.t",
                    help="file with hscore. Default: hscore_file.t")

args = parser.parse_args()

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

# solar options: protein to protein, construct multiblocks, input format m8 (blast tabular)
solar_options = " -a prot2prot -c -f m8"

solar_string = "%s %s %s > %s" % (args.solar_script_path, solar_options, args.input, args.solar_output_file)
print(solar_string)
os.system(solar_string)

with open(args.solar_output_file, "r") as solar_fd:
    with open(args.modified_solar_output_file, "w") as filtered_fd:
        for line in solar_fd:
            line_list = line.strip().split("\t")
            query_id = line_list[0]
            query_len = int(line_list[1])
            seq_id = line_list[5]
            seq_len = int(line_list[6])
            num_of_blocks = int(line_list[9])
            alignment_total_score = line_list[10]
            query_entries = [map(int, entry.split(",")) for entry in line_list[11][:-1].split(";")]
            total_query = sum(map(lambda x: x[1] - x[0] + 1, query_entries))
            sequence_entries = [map(int, entry.split(",")) for entry in line_list[12][:-1].split(";")]
            sequence_query = sum(map(lambda x: x[1] - x[0] + 1, sequence_entries))
            score_entries = line_list[13][:-1].split(";")
            length_ratio = min(float(total_query)/float(query_len), float(sequence_query)/float(seq_len))
            if length_ratio >= 0.33:
                filtered_fd.write("%s\t%i\t%s\t%i\t%s\t%i\t%i\t%f\n" %
                                  (query_id, query_len, seq_id, seq_len, alignment_total_score, total_query,
                                   sequence_query, length_ratio))
"""
get_bit_score_string = "awk -F'\\t' '{printf \"%%s\\t%%s\\t%%s\\n\",$1,$6,$11}' %s > %s" % \
                       (args.solar_output_file, args.bit_score_file)
"""

get_bit_score_string = "awk -F'\\t' '{printf \"%%s\\t%%s\\t%%s\\n\",$1,$3,$5}' %s > %s" % \
                       (args.modified_solar_output_file, args.bit_score_file)

os.system(get_bit_score_string)
self_alignment_score_string = "awk -F'\\t' '{if ($1 == $2) {printf \"%%s\\t%%s\\n\",$1,$3}}' %s > %s" % \
                              (args.bit_score_file, args.self_alignment_score_file)
os.system((self_alignment_score_string))
nodes_dict = {}

with open(args.self_alignment_score_file, "r") as node_fd:
    for line in node_fd:
        name, weight = line.strip().split("\t")
        nodes_dict[name] = int(weight)

with open(args.bit_score_file, "r") as input_fd:
    with open(args.temp_file, "w") as temp_fd:
        #with open(args.output_file, "w") as out_fd:
        for line in input_fd:
            first, second, weight = line.strip().split("\t")
            f, s = sorted([first, second])
            if first == second:
                continue
                #out_fd.write("%s\t%s\t%s\n" % (f, s,  int(100 * float(weight)/ float(max(nodes_dict[first], nodes_dict[second])))))
            temp_fd.write("%s\t%s\t%s\n" % (f, s,  int(100 * float(weight)/float(max(nodes_dict[first], nodes_dict[second])))))

sort_string = "sort -u %s > %s" % (args.temp_file, args.hscore_file)
os.system(sort_string)

hcluster_options = " -w 5 -s 0.33"
hcluster_string = "%s %s %s > %s" % (args.hcluster_sg_path, hcluster_options, args.hscore_file, args.output)

os.system(hcluster_string)
if args.output != "output":
    out_fd.close()
if args.input != "stdin":
    in_fd.close()