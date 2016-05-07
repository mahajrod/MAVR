#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from multiprocessing import Pool


def call_hcluster(input, output, edge_weight, edge_density):
    hcluster_options = " -w %i -s %f" % (edge_weight, edge_density)
    hcluster_string = "%s %s %s > %s" % (args.hcluster_sg_path, hcluster_options, input, output)
    print (hcluster_string)
    os.system(hcluster_string)


def call_hcluster_single_arg(arg_list):
    call_hcluster(args.input, *arg_list)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with hscores. Default: stdin")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", default="gene_clusters",
                    help="Prefix of output file with clustered genes. Default: gene_clusters")
parser.add_argument("-g", "--hcluster_sg_path", action="store", dest="hcluster_sg_path", default="hcluster_sg",
                    help="Path to hcluster_sg binary. Default: hcluster_sg")
parser.add_argument("-e", "--statistics_script_path", action="store", dest="statistics_script_path",
                    default="expansion_pipeline_statistics.py",
                    help="Path to script drawing statistics. Default: expansion_pipeline_statistics.py")
parser.add_argument("-c", "--hscore_file", action="store", dest="hscore_file", default="hscore_file.t",
                    help="file with hscore. Default: hscore_file.t")
parser.add_argument("-w", "--min_weight_of_edge", action="store", dest="min_weight_of_edge", type=int,
                    default=5, help="Starting value of minimum weight of edge. Should be between 0 and 100. Default: 5")
parser.add_argument("-a", "--maximum_weight_of_edge", action="store", dest="max_weight_of_edge", type=int,
                    default=100, help="End value of weight of edge. Should be between 0 and 100. Default: 100")
parser.add_argument("-y", "--step_weight_of_edge", action="store", dest="step_of_weight_of_edge", type=int,
                    default=1, help="Step of weight of edge. Default: 1")

parser.add_argument("-d", "--min_density_of_edge", action="store", dest="min_density_of_edge", type=float,
                    default=0.10, help="Starting value of density of edge. Should be between 0 and 1.00 Default: 0.10")
parser.add_argument("-x", "--max_density_of_edge", action="store", dest="max_density_of_edge", type=float,
                    default=1.00, help="End value of density of edge. Should be between 0 and 1.00 Default: 1.00")
parser.add_argument("-t", "--step_of_density_of_edge", action="store", dest="step_density_of_edge", type=float,
                    default=0.05, help="Step of density of edge. Default: 0.05")
parser.add_argument("-p", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")

parser.add_argument("-n", "--number_of_species", action="store", dest="number_of_species", type=int,
                    help="Number of species in analysis.")
parser.add_argument("-s", "--species_set", action="store", dest="species_set",
                    help="Comma separated set of species.")
parser.add_argument("-l", "--name_last", action="store_false", dest="name_first", default=True,
                    help="Position of name of species in gene_id")


args = parser.parse_args()

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")


print("Running hcluster...")

arguments_list = []
weight = args.min_weight_of_edge
density = args.min_density_of_edge
while weight < args.max_weight_of_edge:
    density = args.min_density_of_edge
    while density < args.max_density_of_edge:
        arguments_list.append(("%s_w_%i_d_%f.t" % (args.output, weight, density), weight, density))
        density += args.step_density_of_edge
    weight += args.step_of_weight_of_edge

processes_pool = Pool(args.threads)
processes_pool.map(call_hcluster_single_arg, arguments_list)

if args.output != "output":
    out_fd.close()
if args.input != "stdin":
    in_fd.close()
