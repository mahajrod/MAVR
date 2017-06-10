#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()
from Bio import SeqIO

from BCBio import GFF

from collections import OrderedDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", default="stdout",
                    help="Output file with distances between genes. Default: stdout")
parser.add_argument("-c", "--coordinates", action="store", dest="coordinates", default="coordinates.bed",
                    help="Output file with coordinates of genes")
parser.add_argument("-g", "--histogram_prefix", action="store", dest="histogram",
                    default="intergenic_distance_distribution",
                    help="Prefix of file with histogram")
parser.add_argument("-f", "--formats", action="store", dest="formats",
                    default=["svg", "png"], type=lambda s: s.split(","),
                    help="Comma-separated list of formats for histogram. Default: svg,png")
parser.add_argument("-b", "--bin_width", action="store", dest="bin_width", type=int,
                    default=1000,
                    help="Width of bin in histogram")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
coordinates_fd = open(args.coordinates, "w")
annotations_dict = SeqIO.to_dict(GFF.parse(open(args.input)))

distance_dict = OrderedDict()
all_distances = []

print(len(annotations_dict))
for record in annotations_dict:
    distances = []
    coordinates = []
    for feature in annotations_dict[record].features:
        #print feature.id
        if feature.type != "gene":
            continue

        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand
        coordinates.append((start, end, feature.id))
        coordinates_fd.write("%s\t%s\t%s\t%s\t%s\n" % (record, start, end, strand, feature.id))

    if len(coordinates) <= 1:
        continue

    coordinates.sort(key=lambda k: k[0])
    for i in range(1, len(coordinates)):
        #if coordinates[i][0] - coordinates[i-1][1] < 0:
            #print record
        distance = coordinates[i][0] - coordinates[i-1][1]
        distances.append(distance)
        out_fd.write("%s\t%i\t%i\t%i\n" % (coordinates[i][2], coordinates[i-1][1], coordinates[i][0], distance))
    all_distances += distances
    distance_dict[record] = np.array(distances)
"""
all_distances.sort()
print len(all_distances)
print all_distances[10150]
print all_distances[9150]
print all_distances[8150]
bins = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000,
        30000,  40000,  50000,  60000,  70000,  80000,
        90000, 100000, 110000, 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000]
print np.mean(np.array(all_distances))
print np.median(np.array(all_distances))
"""
bins = np.arange(-50000, 205000, args.bin_width)
plt.figure(1)
#figure = plt.figure(1, figsize=(30, 30))
subplot = plt.subplot(1, 1, 1)
plt.hist(all_distances, bins=bins)
subplot.tick_params(direction='out')
for filetype in args.formats:
    plt.savefig("%s_bin_width_%i.%s" % (args.histogram, args.bin_width, filetype))


coordinates_fd.close()


"""
sequence_groups_id = SynDict()
sequence_groups_id.read(args.id_file, split_values=True)
#print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)
for group in sequence_groups_id:
    SeqIO.write(record_by_id_generator(sequence_dict, sequence_groups_id[group]),
                "%s%s.%s" % (args.output, group, args.extension), format=args.format)
"""





