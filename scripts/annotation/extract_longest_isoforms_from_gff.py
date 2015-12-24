#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from copy import deepcopy

from BCBio import GFF


def record_with_extracted_annotations_generator(gff_file, mode):
    for record in GFF.parse(open(gff_file)):
        print("Handling region '%s'" % record.id)
        new_record = deepcopy(record)
        new_record.features = []
        for feature in record.features:
            max_len = None
            if feature.type == "gene":
                new_gene_record = deepcopy(feature)
                new_gene_record.sub_features = []
                transcript_record_list = []
                transcript_len_list = []
                CDS_len_list = []
                for subfeature in feature.sub_features:
                    if subfeature.type == "mRNA" or subfeature.type == "transcript":
                        transcript_record_list.append(subfeature)
                        transcript_len_list.append(len(subfeature))
                        CDS_len = 0
                        for subsubfeature in subfeature.sub_features:
                            if subsubfeature.type == "CDS":
                                CDS_len += len(subsubfeature)
                        CDS_len_list.append(CDS_len)
                if mode == "longest_CDS" and CDS_len_list:
                    max_len = max(CDS_len_list)
                    longest = transcript_record_list[CDS_len_list.index(max_len)]
                elif mode == 'longest_transcript' and transcript_len_list:
                    max_len = max(transcript_len_list)
                    longest = transcript_record_list[transcript_len_list.index(max_len)]
                if max_len:
                    new_gene_record.sub_features.append(longest)
                    if max_len > 0:
                        new_record.features.append(new_gene_record)
        if new_record.features:
            yield new_record

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Output file with longest transcripts")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="longest_CDS",
                    help="Operational mode. Possible variants: 'longest_CDS', 'longest_transcript'. "
                         "Default - longest_CDS")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

GFF.write(record_with_extracted_annotations_generator(args.input_gff, args.mode), out_fd)

out_fd.close()
