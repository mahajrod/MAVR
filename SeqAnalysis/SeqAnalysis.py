#!/usr/bin/env python2

import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
pyplot.ioff()
from random import randint

from Bio import SeqIO

from Bio.SeqFeature import SeqFeature, FeatureLocation


def draw_species_counts_distribution(number_of_species_dict,
                                     distribution_file,
                                     filter_min_value=None,
                                     filter_max_value=None):

    figure = pyplot.figure(dpi=400, figsize=(30, 20))
    axes = figure.add_subplot(1, 1, 1)
    values = []
    print("Drawning distribution of species counts...")
    if filter_min_value and filter_max_value:
        for value in number_of_species_dict.values():
            if value >= filter_min_value  and value <= filter_max_value:
                values.append(value)
        print("Counting only species with sequences in range [%i, %i]" % (filter_min_value, filter_max_value))
    elif filter_min_value:
        for value in number_of_species_dict.values():
            if value >= filter_min_value:
                values.append(value)
        print("Counting only species with % or more sequences" % filter_min_value)
    elif filter_max_value:
        for value in number_of_species_dict.values():
            if value <= filter_max_value:
                values.append(value)
        print("Counting only species with % or less sequences" % filter_max_value)
    else:
        values = list(number_of_species_dict.values())
        print("No filter was set for number of sequences per species")
    max_value = max(values)
    print(values, max_value)
    n, bins, patches = axes.hist(values, max_value, facecolor='green')
    axes.set_xlabel('Number of sequences per species')
    axes.set_ylabel('Number of species')
    axes.set_title('Distribution of sequences per species')
    print("Totaly %i species" % len(values))
    pyplot.savefig(distribution_file)
    axes.set_yscale('log')
    pyplot.savefig('log_' + distribution_file)

def write_sequences_coordinates_to_file(seq_coordinates_tuple_list, output_filename):
    fd = open(output_filename, "w")
    fd.write("#length\tstart\tend\n")
    for coord_tuple in seq_coordinates_tuple_list:
        fd.write("%i\t%i\t%i\n" % (coord_tuple[1] - coord_tuple[0] + 1, coord_tuple[0], coord_tuple[1]))
    fd.close()
    return 1

def replace_ambigious_nucleotides(input_file, output_file, index_file="index.idx", format="fasta"):
    record_dict = SeqIO.index_db(index_file, [input_file], format)
    corrected_dict = {}
    ambigious_reg_exp = re.compile("[^AGTCN]+", re.IGNORECASE)
    ambigious_dict = {}
    for record in record_dict:
        ambigious_dict[record] = []
        corrected_dict[record] = record_dict[record]
        #print self.reference_genome[entry].seq
        ambigious = ambigious_reg_exp.finditer(str(record_dict[record].seq))  # iterator with
        for match in ambigious:
            #print(match)
            ambigious_dict[record].append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                         type="gap", strand=None))
