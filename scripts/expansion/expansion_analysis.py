#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with blast alignment. Default: stdin")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", default="gene_clusters",
                    help="Prefix of output file with clustered genes. Default: gene_clusters")
parser.add_argument("-a", "--solar_script_path", action="store", dest="solar_script_path", default="solar.pl",
                    help="Path to solar script. Default: solar.pl")
parser.add_argument("-g", "--hcluster_sg_path", action="store", dest="hcluster_sg_path", default="hcluster_sg",
                    help="Path to hcluster_sg binary. Default: hcluster_sg")
parser.add_argument("-l", "--solar_output_file", action="store", dest="solar_output_file", default="solar_output.t",
                    help="File with solar output. Default: solar_output.t")
parser.add_argument("-b", "--bit_score_file", action="store", dest="bit_score_file", default="bit_score.t",
                    help="File with solar bitscores. Default: bit_score.t")
parser.add_argument("-e", "--self_alignment_score_file", action="store", dest="self_alignment_score_file", default="self_alignment_score.t",
                    help="File with self alignment scores. Default: self_alignment_score.t")
parser.add_argument("-t", "--temp_file", action="store", dest="temp_file", default="temp_file.t",
                    help="Temp file. Default: temp_file.t")
parser.add_argument("-c", "--hscore_file", action="store", dest="hscore_file", default="hscore_file.t",
                    help="file with hscore. Default: hscore_file.t")
parser.add_argument("-s", "--size_of_common_fragment", action="store", dest="size_of_common_fragment", type=float,
                    default=0.33, help="Minimum threshold of size of fragment mapped to both compared genes. Should be between 0 and 1.00 Default: 0.33")
parser.add_argument("-w", "--minimum_weight_of_edge", action="store", dest="minimum_weight_of_edge", type=int,
                    default=5, help="Minimum weight of edge. Should be between 0 and 100. Default: 5")
parser.add_argument("-d", "--minimum_density_of_edge", action="store", dest="minimum_density_of_edge", type=float,
                    default=0.33, help="Minimum density of edge. Should be between 0 and 1.00 Default: 0.33")


args = parser.parse_args()

if args.size_of_common_fragment > 1:
    raise ValueError("Size of common fragment should between 0 and 1.00")
elif (not isinstance(args.minimum_weight_of_edge, int)) or (args.minimum_density_of_edge > 100):
    raise ValueError("Minimum weight of edge should be integer and less then 100")
elif args.minimum_density_of_edge > 1:
    raise ValueError("Minimum density of edge should between 0 and 1.00")

# solar options: protein to protein, nput format m8 (blast tabular)
solar_options = " -a prot2prot -f m8 -z"

solar_string = "%s %s %s > %s" % (args.solar_script_path, solar_options, args.input, args.solar_output_file)
print("Clustering segments...\n\t" + solar_string)
os.system(solar_string)

print("Extracting bit scores...")
nodes_dict = {}
length_dict = {}
with open(args.solar_output_file, "r") as solar_fd:
    with open(args.bit_score_file, "w") as filtered_fd:
        with open(args.self_alignment_score_file, "w") as self_aligned_fd:
            for line in solar_fd:
                line_list = line.strip().split("\t")
                query_id = line_list[0]
                query_len = int(line_list[1])
                seq_id = line_list[5]
                #seq_len = int(line_list[6])
                #num_of_blocks = int(line_list[9])
                alignment_total_score = line_list[10]
                if query_id == seq_id:
                    self_aligned_fd.write("%s\t%i\t%s\n" % (query_id, query_len, alignment_total_score))
                    nodes_dict[query_id] = int(alignment_total_score)
                    length_dict[query_id] = query_len
                else:
                    query_entries = [map(int, entry.split(",")) for entry in line_list[11][:-1].split(";")]
                    total_query = sum(map(lambda x: x[1] - x[0] + 1, query_entries)) # total length of alignned segments of query
                    sequence_entries = [map(int, entry.split(",")) for entry in line_list[12][:-1].split(";")]
                    sequence_query = sum(map(lambda x: x[1] - x[0] + 1, sequence_entries)) # total length of sequence segments corresponding to alignment
                    alignment_length = min((total_query, sequence_query))
                    f, s = sorted([query_id, seq_id])
                    filtered_fd.write("%s\t%s\t%s\t%i\n" %
                                      (f, s, alignment_total_score, alignment_length))

sort_string = "sort %s > %s" % (args.bit_score_file, args.temp_file)
print("Sorting...\n\t%s" % sort_string)
os.system(sort_string)

print("Calculating hscore...")
with open(args.temp_file, "r") as input_fd:
    with open(args.hscore_file, "w") as hscore_fd:
        #with open(args.output_file, "w") as out_fd:
        for line in input_fd:
            fl_first, fl_second, fl_weight, fl_alignment_len = line.strip().split("\t")
            fl_alignment_len = float(fl_alignment_len)
            try:  # check for last string in file (Stop Iteration is raised if end of file was reached on .next() method)
                sl_first, sl_second, sl_weight, sl_alignment_len = input_fd.next().strip().split("\t")
            except StopIteration:
                break
            sl_alignment_len = float(sl_alignment_len)
            while (fl_first != sl_first) or (fl_second != sl_second):
                fl_first, fl_second, fl_weight, fl_alignment_len = sl_first, sl_second, sl_weight, sl_alignment_len
                sl_first, sl_second, sl_weight, sl_alignment_len = input_fd.next().strip().split("\t")
                sl_alignment_len = float(sl_alignment_len)
            if (fl_first not in length_dict) or (fl_second not in length_dict):
                continue
            length_ratio = min([fl_alignment_len, sl_alignment_len]) / max([length_dict[fl_first],
                                                                            length_dict[fl_second]])
            """
            print(fl_first, fl_second, fl_weight, fl_alignment_len)
            print(sl_first, sl_second, sl_weight, sl_alignment_len)
            print(length_dict[fl_first], length_dict[fl_second])
            print(length_ratio)
            """
            if length_ratio > args.size_of_common_fragment:
                fl_score = int(100 * float(fl_weight)/float(max([nodes_dict[fl_first], nodes_dict[fl_second]])))
                sl_score = int(100 * float(sl_weight)/float(max([nodes_dict[sl_first], nodes_dict[sl_second]])))
                hscore = min([fl_score, sl_score])
                hscore_fd.write("%s\t%s\t%s\n" % (fl_first, fl_second, hscore))

hcluster_options = " -w %i -s %f" % (args.minimum_weight_of_edge, args.minimum_density_of_edge)
hcluster_string = "%s %s %s > %s" % (args.hcluster_sg_path, hcluster_options, args.hscore_file,
                                     "%s_w_%i_d_%f.t" % (args.output, args.minimum_weight_of_edge,
                                                         args.minimum_density_of_edge))
print("Clustering...\n\t%s" % hcluster_string)
os.system(hcluster_string)

os.remove(args.temp_file)