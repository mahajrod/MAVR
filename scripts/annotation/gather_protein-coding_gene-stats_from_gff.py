#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from copy import deepcopy
from collections import OrderedDict
from BCBio import GFF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations of protein-coding genes")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with statistics")

args = parser.parse_args()


out_fd = open(args.output_file, "w")
in_fd = open(args.input_gff, "r")

out_fd.write("#gene\tN_of_transcripts\tTranscripts_ids\tLength_of_transcripts\tN_of_exons\tLength_of_exons\tLength_of_CDS\n")

for record in GFF.parse(in_fd):
    for feature in record.features:
        mRNA_exon_dict = OrderedDict()
        mRNA_CDS_dict = OrderedDict()
        mRNA_intron_dict = OrderedDict()
        mRNA_length_dict = OrderedDict()

        if feature.type != "gene":
            continue
        print("Handling %s" % feature.id)
        for subfeature in feature.sub_features:
            if subfeature.type != "mRNA" or subfeature.type != "transcript":
                continue
            mRNA_length_dict[subfeature.id] = 0
            mRNA_exon_dict[subfeature.id] = []
            mRNA_intron_dict[subfeature.id] = []
            mRNA_CDS_dict[subfeature.id] = []

            for subsubfeature in subfeature.sub_features:
                if subsubfeature.type == "exon":
                    mRNA_exon_dict[subfeature.id].append(len(subsubfeature))
                    mRNA_length_dict[subfeature.id] += len(subsubfeature)
                elif subsubfeature.type == "CDS":
                    mRNA_CDS_dict[subfeature.id].append(len(subsubfeature))

        #mRNA_len_str = ""
        mRNA_len_str = ";".join(map(str, mRNA_length_dict.values()))
        mRNA_N_of_exons = ";".join([str(len(mRNA_exon_dict[mRNA]))for mRNA in mRNA_exon_dict])
        mRNA_length_of_exons = ";".join([",".join(map(str, mRNA_exon_dict[mRNA])) for mRNA in mRNA_exon_dict])
        mRNA_length_of_CDS = ";".join([",".join(map(str, mRNA_CDS_dict[mRNA])) for mRNA in mRNA_CDS_dict])
        out_fd.write("%s\t%i\t%s\t" % (feature.id, len(mRNA_exon_dict), ";".join(mRNA_exon_dict.keys())))
        out_fd.write("%s\t%s\t%s\t%s\n" % (mRNA_len_str, mRNA_N_of_exons, mRNA_length_of_exons, mRNA_length_of_CDS))

in_fd.close()
out_fd.close()
