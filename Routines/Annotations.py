__author__ = 'mahajrod'

from copy import deepcopy

from Bio import SearchIO
from BCBio import GFF

from CustomCollections.GeneralCollections import IdSet, SynDict, IdList


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

    @staticmethod
    def replace_region_names_in_gff(input_gff, synonyms_file, output_gff):
        syn_dict = SynDict()
        syn_dict.read(synonyms_file, comments_prefix="#")
        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                    else:
                        line_list = line.split("\t")
                        if line_list[0] in syn_dict:
                            line_list[0] = syn_dict[line_list[0]]
                            out_fd.write("\t".join(line_list))
                        else:
                            out_fd.write(line)

    @staticmethod
    def get_transcript_to_pep_accordance_from_gtf(gtf_file, output_file, comment_symbol="#"):
        """
        Tested on gtf files from Ensembl relealese 70
        """
        accordance_dict = SynDict()
        with open(gtf_file, "r") as gtf_fd:
            for line in gtf_fd:
                if line[0] == comment_symbol:
                    continue
                tmp_list = line.strip().split("\t")
                tmp_list = tmp_list[-1].split(";")
                protein_id = None
                transcript_id = None
                #print tmp_list
                for entry in tmp_list:
                    tmp_entry = entry.split()

                    if len(tmp_entry) != 2:
                        continue
                    if tmp_entry[0] == "transcript_id":
                        #print "tttt"
                        transcript_id = tmp_entry[1][1:-1]  # remove quotes
                    elif tmp_entry[0] == "protein_id":
                        #print "ppppp"
                        protein_id = tmp_entry[1][1:-1]

                if (transcript_id is not None) and (protein_id is not None):
                    if transcript_id in accordance_dict:
                        accordance_dict[transcript_id].add(protein_id)
                    else:
                        accordance_dict[transcript_id] = {protein_id}
        accordance_dict.write(output_file, splited_values=True)

    @staticmethod
    def fix_gff_coordinates_order(input_gff, output_gff):
        fixed_lines_number = 0
        malformed_lines_number = 0
        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    tmp = line.split("\t")

                    try:
                        start = int(tmp[3])
                        end = int(tmp[4])
                    except ValueError:
                        malformed_lines_number += 1
                        continue

                    if start <= end:
                        out_fd.write(line)
                        continue

                    tmp[3], tmp[4] = tmp[4], tmp[3]
                    out_fd.write("\t".join(tmp))

                    fixed_lines_number += 1

        print("Fixed lines: %i" % fixed_lines_number)
        print("Malformed lines(removed): %i" % malformed_lines_number)

    @staticmethod
    def fix_absent_feature_type_field(input_gff, output_gff, feature_type):
        fixed_lines_number = 0
        malformed_lines_number = 0
        malformed_lines_list = []
        line_number = 0
        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    line_number += 1
                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    tmp = line.split("\t")

                    try:
                        start = int(tmp[3])
                        end = int(tmp[4])

                    except ValueError:
                        malformed_lines_number += 1
                        malformed_lines_list.append(line_number)
                        tmp = tmp[:2] + [feature_type] + tmp[2:]
                        out_fd.write("\t".join(tmp))
                        continue

                    out_fd.write(line)

    @staticmethod
    def add_alias_to_feature(input_gff, output_gff, syn_file, feature_type_list=(), name_field_list=("ID",),
                             alias_field="Alias", key_column=0, value_column=1):
        syn_dict = SynDict(filename=syn_file, comments_prefix="#", split_values=True,
                           values_separator=",", key_index=key_column, value_index=value_column)
        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:

                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    tmp = line.strip().split("\t")
                    if feature_type_list:  # skip feature types not present in feature type list
                        if tmp[2] not in feature_type_list:
                            out_fd.write(line)
                            continue
                    description_list = tmp[8].split(";")

                    feature_aliases = []
                    for description_entry in description_list:
                        entry_list = description_entry.split("=")
                        if entry_list[0] not in name_field_list:
                            continue
                        feature_aliases += entry_list[1].split(",")

                    if not feature_aliases:
                        print("No ID was found for this line: \n%s" % line)

                    feature_synonyms_list = []
                    syn_lines = 0
                    for name in feature_aliases:
                        if name in syn_dict:
                            feature_synonyms_list += syn_dict[name]
                            syn_lines += 1
                    #print syn_dict
                    if syn_lines > 1:
                        print "Warning!!!. For gene with names %s were found more 1 synonym group." % ",".join(feature_aliases)

                    for i in range(0, len(description_list)):
                        entry_list = description_list[i].split("=")
                        if entry_list[0] != alias_field:
                            continue

                        alias_list = set(entry_list[1].split(",") + feature_synonyms_list)
                        description_list[i] = "%s=%s" % (entry_list[0], ",".join(alias_list))
                    else:
                        description_list.append("%s=%s" % (alias_field, ",".join(feature_synonyms_list)))

                    out_fd.write("\t".join(tmp[:8] + [";".join(description_list)]) + "\n")

    @staticmethod
    def count_total_feature_length_from_gff(input_gff, output_prefix, features_to_count=None):
        len_file = "%s.len" % output_prefix
        stat_file = "%s.stat" % output_prefix

        total_feature_length = 0
        feature_number = 0

        with open(input_gff, "r") as in_fd:
            with open(len_file, "w") as len_fd:
                for line in in_fd:
                    if line[0] == "#":
                        continue
                    tmp = line.split("\t")
                    feature = tmp[2]
                    if features_to_count is not None:
                        if feature not in features_to_count:
                            continue

                    start = int(tmp[3])
                    end = int(tmp[4])

                    feature_number += 1
                    feature_length = end - start + 1

                    len_fd.write("%i\n" % feature_length)

                    total_feature_length += feature_length

        stat_string = "Features\t%s\n" % (",".join(features_to_count) if features_to_count else "all")
        stat_string += "Number of features\t%i\n" % feature_number
        stat_string += "Total length\t%i\n" % total_feature_length

        print(stat_string)
        with open(stat_file, "w") as stat_fd:
            stat_fd.write(stat_string)

    @staticmethod
    def get_feature_to_parent_correspondence_from_gff(input_gff, output, feature_list=("mRNA",),
                                                      id_entry="ID", parental_id_entry="Parent"):
        with open(input_gff, "r") as in_fd:
            with open(output, "w") as out_fd:
                for line in in_fd:
                    if line[0] == "#":
                        continue
                    tmp = line.strip().split("\t")
                    if tmp[2] not in feature_list:
                        continue
                    tmp = map(lambda s: s.split("="), tmp[-1].split(";"))

                    feature_id = None
                    parental_id = None

                    for entry, value in tmp:
                        if entry == id_entry:
                            feature_id = value
                        elif entry == parental_id_entry:
                            parental_id = value
                    if (parental_id is None) or (feature_id is None):
                        print ("Absent feature id or parental id in following line:")
                        print (line)
                        continue
                    out_fd.write("%s\t%s\n" % (parental_id, feature_id))

    @staticmethod
    def add_length_to_accordance_file(accordance_file, length_file, output_prefix):

        accordance_dict = SynDict(filename=accordance_file)
        length_dict = SynDict(length_file, expression=int)
        longest_list = IdList()

        all_output_file = "%s.all.correspondce" % output_prefix
        longest_output_file = "%s.longest.correspondence" % output_prefix
        longest_id_file = "%s.longest.ids" % output_prefix

        current_gene = None
        current_transcript = None
        current_length = 0

        with open(all_output_file, "w") as all_out_fd:
            with open(longest_output_file, "w") as longest_out_fd:
                for gene in accordance_dict:
                    all_out_fd.write("%s\t%s\t%i\n" % (gene, accordance_dict[gene], length_dict[accordance_dict[gene]]))

                    if current_gene and (gene != current_gene):
                        longest_out_fd.write("%s\t%s\t%i\n" % (current_gene, current_transcript, current_length))
                        longest_list.append(current_transcript)
                        current_gene = gene
                        current_transcript = accordance_dict[gene]
                        current_length = length_dict[accordance_dict[gene]]
                        continue

                    if length_dict[accordance_dict[gene]] > current_length:
                        current_gene = gene
                        current_transcript = accordance_dict[gene]
                        current_length = length_dict[accordance_dict[gene]]











