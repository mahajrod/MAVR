#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from BCBio import GFF

from CustomCollections.GeneralCollections import IdList


def record_with_extracted_annotations_generator(gff_file):
    for record in GFF.parse(open(gff_file)):
        print("Extracting annotations from %s" % record.id)
        new_record = record
        record.features = []
        for feature in record.features:
            if args.types:
                if feature.id in annotation_ids:
                    new_record.features.append(feature)
        yield new_record

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with extracted_annotations")
parser.add_argument("-d", "--ids_file", action="store", dest="ids_file",
                    help="File with ids of annotations to extract")

args = parser.parse_args()

annotation_ids = IdList()
annotation_ids.read(args.ids_file, comments_prefix="#")

GFF.write(record_with_extracted_annotations_generator(args.input_gff), open(args.output_file, "w"))
