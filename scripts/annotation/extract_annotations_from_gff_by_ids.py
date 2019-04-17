#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from copy import deepcopy
from BCBio import GFF
from RouToolPa.Collections.General import IdList



def record_with_extracted_annotations_generator(gff_file, white_list_of_annotation_types):
    for record in GFF.parse(open(gff_file)):
        #print("Extracting annotations from %s" % record.id)
        new_record = deepcopy(record)
        new_record.features = []
        #print record.features
        for feature in record.features:

            #print ("%s\t%s" % (record.id, feature.id))

            if (feature.id in annotation_ids) and (feature.type in white_list_of_annotation_types):
                new_record.features.append(feature)
        if len(new_record.features) > 0:
            yield new_record

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with extracted_annotations")
parser.add_argument("-d", "--ids_file", action="store", dest="ids_file",
                    help="File with ids of annotations to extract")
parser.add_argument("-t", "--annotation_types", action="store", dest="annotation_types", default=["gene"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of annotation types to extract")

args = parser.parse_args()

annotation_ids = IdList()
annotation_ids.read(args.ids_file, comments_prefix="#")
#print args.annotation_types
out_fd = open(args.output_file, "w")

GFF.write(record_with_extracted_annotations_generator(args.input_gff, args.annotation_types), out_fd)

out_fd.close()
