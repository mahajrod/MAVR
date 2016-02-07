__author__ = 'mahajrod'

from copy import deepcopy

from Bio import SearchIO
from BCBio import GFF

from CustomCollections.GeneralCollections import IdSet


class AnnotationsRoutines:
    def __init__(self):
        pass

    @staticmethod
    def record_with_extracted_annotations_generator(gff_file,
                                                    annotation_ids,
                                                    white_list_of_annotation_types):
        for record in GFF.parse(open(gff_file)):
            new_record = deepcopy(record)
            new_record.features = []
            for feature in record.features:
                if (feature.id in annotation_ids) and (feature.type in white_list_of_annotation_types):
                    new_record.features.append(feature)
            if len(new_record.features) > 0:
                yield new_record

    def extract_annotation_from_gff(self, input_gff, annotation_ids,
                                    white_list_of_annotation_types, output_gff):

        GFF.write(self.record_with_extracted_annotations_generator(input_gff, annotation_ids,
                                                                   white_list_of_annotation_types),
                  open(output_gff, "w"))

    @staticmethod
    def record_with_extracted_transcripts_generator(gff_file,
                                                    transcript_ids):
        for record in GFF.parse(open(gff_file)):
            new_record = deepcopy(record)
            new_record.features = []
            for feature in record.features:
                if (feature.type == "mRNA" or feature.type == "transcript") and (feature.id in transcript_ids):
                    new_record.features.append(feature)
                elif feature.type == "gene":
                    new_feature = deepcopy(feature)
                    new_feature.sub_features = []

                    for subfeature in feature.sub_features:
                        if (subfeature.type == "mRNA" or subfeature.type == "transcript") and (subfeature.id in transcript_ids):
                            new_feature.sub_features.append(subfeature)
                    if len(new_feature.sub_features) > 0:
                        new_record.features.append(new_feature)
            if len(new_record.features) > 0:
                yield new_record

    def extract_transcripts_by_ids(self, input_gff, transcript_id_file, output_gff):
        transcript_ids = IdSet()
        transcript_ids.read(transcript_id_file, header=False)
        GFF.write(self.record_with_extracted_transcripts_generator(input_gff, transcript_ids),
                  open(output_gff, "w"))
