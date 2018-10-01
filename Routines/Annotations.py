__author__ = 'mahajrod'
import os
from copy import deepcopy

from collections import OrderedDict
from Bio import SearchIO
from BCBio import GFF

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from Routines.Sequence import SequenceRoutines
from CustomCollections.GeneralCollections import IdSet, SynDict, IdList


class AnnotationsRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)

        self.gff_scaffold_column = 0
        self.gff_source_column = 1
        self.gff_featuretype_column = 2
        self.gff_start_column = 3
        self.gff_end_column = 4
        self.gff_score_column = 5
        self.gff_strand_column = 6
        self.gff_phase_column = 7
        self.gff_attribute_column = 8


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
                        print("Warning!!!. For gene with names %s were found more 1 synonym group." % ",".join(feature_aliases))

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

        accordance_dict = SynDict(filename=accordance_file, allow_repeats_of_key=True)
        length_dict = SynDict(filename=length_file, expression=int)
        print length_dict
        longest_list = IdList()

        all_output_file = "%s.all.correspondence" % output_prefix
        longest_output_file = "%s.longest.correspondence" % output_prefix
        longest_id_file = "%s.longest.ids" % output_prefix

        with open(all_output_file, "w") as all_out_fd:
            with open(longest_output_file, "w") as longest_out_fd:
                for gene in accordance_dict:
                    current_transcript = None
                    current_length = 0
                    for transcript in accordance_dict[gene]:
                        if length_dict[transcript] > current_length:
                            current_transcript = transcript
                            current_length = length_dict[transcript]
                        all_out_fd.write("%s\t%s\t%i\n" % (gene, transcript, length_dict[transcript]))

                    longest_out_fd.write("%s\t%s\t%i\n" % (gene, current_transcript, current_length))
                    longest_list.append(current_transcript)
        longest_list.write(longest_id_file)

    @staticmethod
    def parse_gff_annotation_string_to_dict(gff_annotation_string):

        annotation_dict = OrderedDict()
        for entry in gff_annotation_string.split(";"):
            key, value = entry.split("=")
            annotation_dict[key] = value.split(",")

        return annotation_dict

    def extract_sequences_by_gff(self, input_file, gff_file, output_file, type_list=("gene",), parsing_mode="parse",
                                 tmp_index_file="temp.idx", format="fasta"):

        annotations_dict = SeqIO.to_dict(GFF.parse(open(gff_file)))
        #print annotations_dict
        #print("Parsing %s..." % args.input)
        sequence_dict = self.parse_seq_file(input_file, parsing_mode, format, index_file=tmp_index_file)  # SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)

        SeqIO.write(self.record_generator(annotations_dict, sequence_dict, type_list), output_file,
                    format=format)

        if parsing_mode == "index_db":
            os.remove(tmp_index_file)

    def extract_gff_records_by_description_value(self, input_gff, output_gff, field_id_list, value_list,
                                                retain_comments=False):
        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    if line[0] == "#":
                        if retain_comments:
                            out_fd.write(line)
                        continue
                    description_dict = self.get_description_dict_from_gff_string(line, split_values=True, value_separator=",")
                    found = False
                    for field_id in field_id_list:
                        if field_id in description_dict:
                            #print "aaaaa"
                            for value in description_dict[field_id]:
                                if value in value_list:
                                    out_fd.write(line)
                                    found = True
                                    break
                            if found:
                                break

    @staticmethod
    def get_description_dict_from_gff_string(gff_string, split_values=False, value_separator=","):
        if gff_string[0] == "#":
            return None

        description_list = gff_string.strip().split("\t")[-1].split(";")
        description_dict = OrderedDict()

        for entry in description_list:
            key, value = entry.split("=")
            description_dict[key] = value.split(value_separator) if split_values else value

        return description_dict

    def filter_gff_by_description(self, input_gff, output_gff, filtered_out_gff, expression):

        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                with open(filtered_out_gff, "w") as filtered_out_fd:
                    for line in in_fd:
                        if line[0] == "#":
                            out_fd.write(line)
                            continue

                        #print line
                        #print expression(self.get_description_dict_from_gff_string(line))

                        if expression(self.get_description_dict_from_gff_string(line)):
                            out_fd.write(line)
                        else:
                            filtered_out_fd.write(line)

    def add_flanks_to_gff_record(self, input_gff, output_prefix, left_flank_len, right_flank_len, fasta_file,
                                 coords_description_entry="core_seq_coords", id_description_entry="ID"):
        sequence_length_dict = self.get_lengths_from_seq_file(fasta_file)
        shorter_flanks_dict = SynDict()

        output_gff = "%s.gff" % output_prefix
        short_flanks_file = "%s.short_flanks.dat" % output_prefix

        with open(input_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    line_list = line.strip().split("\t")
                    scaffold = line_list[0]
                    start = int(line_list[3])
                    end = int(line_list[4])

                    record_id = OrderedDict(map(lambda s: s.split("="), line_list[8].split(";")))[id_description_entry]

                    line_list[8] += ";%s=%i,%i" % (coords_description_entry, start, end)

                    if line_list[6] == "-":
                        if start - right_flank_len > 0:
                            line_list[3] = str(start - right_flank_len)
                            right_flank_length = right_flank_len
                        else:
                            right_flank_length = start - 1
                            line_list[3] = "1"

                        if end + left_flank_len <= sequence_length_dict[line_list[0]]:
                            line_list[4] = str(end + left_flank_len)
                            left_flank_length = left_flank_len
                        else:
                            left_flank_length = sequence_length_dict[line_list[0]] - end
                            line_list[4] = sequence_length_dict[line_list[0]]
                    else:
                        if start - left_flank_len > 0:
                            line_list[3] = str(start - left_flank_len)
                            left_flank_length = left_flank_len
                        else:
                            left_flank_length = start - 1
                            line_list[3] = "1"

                        if end + right_flank_len <= sequence_length_dict[line_list[0]]:
                            line_list[4] = str(end + right_flank_len)
                            right_flank_length = right_flank_len
                        else:
                            right_flank_length = sequence_length_dict[line_list[0]] - end
                            line_list[4] = str(sequence_length_dict[line_list[0]])

                    if (left_flank_length < left_flank_len) or (right_flank_length < right_flank_len):
                        print("%s: Short flank" % record_id)
                        shorter_flanks_dict[record_id] = "%i,%i" % (left_flank_length, right_flank_length)
                    line_list[8] += ";%s_relative=%i,%i\n" % (coords_description_entry,
                                                              1 + (right_flank_length if line_list[6] == "-" else left_flank_length),
                                                              end - start + 1 + (right_flank_length if line_list[6] == "-" else left_flank_length))
                    """
                    print line
                    print line_list
                    for element in line_list:
                        print element
                        print type(element)
                    """
                    out_fd.write("\t".join(line_list))

        shorter_flanks_dict.write(short_flanks_file)

    def convert_gff_to_bed(self, input_gff, output_bed, id_entry="ID", feature_type_list=[]):
        with open(input_gff, "r") as gff_fd:
            with open(output_bed, "w") as bed_fd:
                for line in gff_fd:
                    if line[0] == "#":
                        continue
                    tmp_list = line.strip().split("\t")
                    if feature_type_list:
                        if tmp_list[2] not in feature_type_list:
                            continue
                    description_dict = self.parse_gff_annotation_string_to_dict(tmp_list[-1])

                    bed_fd.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (tmp_list[0],
                                                               tmp_list[3],
                                                               tmp_list[4],
                                                               description_dict[id_entry],
                                                               tmp_list[5],
                                                               tmp_list[6]))

    @staticmethod
    def convert_gff_to_simple_bed(input_gff, output_bed, feature_type_list=[], scaffold_id_file=None):
        if scaffold_id_file:
            scaffolds_id_list = IdList(filename=scaffold_id_file)
        with open(input_gff, "r") as gff_fd:
            with open(output_bed, "w") as bed_fd:
                for line in gff_fd:
                    if line[0] == "#":
                        continue
                    tmp_list = line.strip().split("\t")
                    if scaffold_id_file:
                        if tmp_list[0] not in scaffolds_id_list:
                            continue
                    if feature_type_list:
                        if tmp_list[2] not in feature_type_list:
                            continue
                    bed_fd.write("%s\t%s\t%s\n" % (tmp_list[0], tmp_list[3], tmp_list[4]))

    def get_id_based_dict_from_gff(self, input_gff, id_entry="ID"):
        id_based_dict = OrderedDict()

        with open(input_gff, "r") as gff_fd:
            for line in gff_fd:
                if line[0] == "#":
                    continue
                tmp_list = line.strip().split("\t")
                description_dict = self.parse_gff_annotation_string_to_dict(tmp_list[-1])

                #print description_dict[id_entry]
                id_based_dict[description_dict[id_entry][0]] = [tmp_list[0], tmp_list[3], tmp_list[4], description_dict]

        return id_based_dict

    @staticmethod
    def get_scaffold_ids_from_gff(gff_file, out_file=None):
        scaffold_id_set = IdSet()

        with open(gff_file, "r") as gff_fd:
            for line in gff_fd:
                if line[0] == "#":
                    continue
                scaffold_id = line.split("\t")[0]
                scaffold_id_set.add(scaffold_id)

        if out_file:
            scaffold_id_set.write(out_file)

        return scaffold_id_set

    @staticmethod
    def count_per_scaffold_feature_number(gff_file, out_file=None, feature_type_list=[]):
        feature_count_dict = SynDict()

        if feature_type_list:
            def check_feature_type(feature_type):
                return feature_type in feature_type_list
        else:
            def check_feature_type(feature_type):
                return True

        with open(gff_file, "r") as gff_fd:
            for line in gff_fd:
                if line[0] == "#":
                    continue
                line_list = line.split("\t")
                if check_feature_type(line_list[2]):
                    if line_list[0] in feature_count_dict:
                        feature_count_dict[line_list[0]] += 1
                    else:
                        feature_count_dict[line_list[0]] = 1

        if out_file:
            feature_count_dict.write(out_file)

        return feature_count_dict

    @staticmethod
    def check_chunks(chunk_dir, number_of_chunks, minimum_chunk_size, separator="_",
                     chunk_filename_suffix=None, chunk_filename_prefix=None):
        chunk_files = os.listdir(chunk_dir)
        filtered_chunk_files = []
        too_small_chunks = []
        chunk_numbers = []

        suffix_length = len(chunk_filename_suffix) if chunk_filename_suffix else None
        prefix_length = len(chunk_filename_prefix) if chunk_filename_prefix else None

        for filename in chunk_files:
            filename_length = len(filename)
            full_filename = "%s/%s" % (chunk_dir, filename)

            if os.path.isdir(full_filename):
                continue

            if chunk_filename_suffix and chunk_filename_prefix:
                if (filename_length <= suffix_length) or (filename_length <= prefix_length):
                    continue
                if (filename[0:prefix_length] == chunk_filename_prefix) and (filename[-suffix_length:] == chunk_filename_suffix):
                    filtered_chunk_files.append(filename)

            elif chunk_filename_prefix:
                if filename_length <= prefix_length:
                    continue
                if filename[0:prefix_length] == chunk_filename_prefix:
                    filtered_chunk_files.append(filename)

            elif chunk_filename_suffix:
                if filename_length <= suffix_length:
                    continue
                if filename[-suffix_length:] == chunk_filename_suffix:
                    filtered_chunk_files.append(filename)
            else:
                filtered_chunk_files.append(filename)

        for filename in filtered_chunk_files:
            full_filename = "%s/%s" % (chunk_dir, filename)
            chunk_number = int(filename.strip(chunk_filename_suffix).split(separator)[-1]) if chunk_filename_suffix else int(filename.split(separator)[-1])
            chunk_numbers.append(chunk_number)

            if os.path.getsize(full_filename) < minimum_chunk_size:
                too_small_chunks.append(chunk_number)

        chunk_numbers.sort()
        too_small_chunks.sort()

        if chunk_numbers[-1] > number_of_chunks:
            print("Largest number of present chunks(%i) is larger than expected(%i)" % (chunk_numbers[-1],
                                                                                        number_of_chunks))

        absent_chunks = []
        for i in range(1, number_of_chunks + 1):
            if i not in chunk_numbers:
                absent_chunks.append(i)

        if absent_chunks:
            print("Absent chunks(total %i): " % len(absent_chunks))
            for k in absent_chunks:
                print(k)
        else:
            print("No absent chunks")

        if too_small_chunks:
            print("Too small chunks (total %i): " % len(too_small_chunks))
            for k in too_small_chunks:
                print(k)
        else:
            print("No too small chunks")

    def merge_overlapping_feature_in_simple_format(self,
                                                   input_file_file_list, scaffold_id_column,
                                                   feature_start_column, feature_end_column,
                                                   output_file=None, output_separator="\t",
                                                   comments_prefix="#", input_separator="\t",
                                                   coordinates_type="1-based", return_seqfeature_dict=False,
                                                   feature_type=None):

        file_list = [input_file_file_list] if isinstance(input_file_file_list, str) else input_file_file_list

        record_dict_list = []

        for filename in file_list:
            record_dict_list.append(OrderedDict())
            for line_list in self.file_line_as_list_generator(filename,
                                                              comments_prefix=comments_prefix,
                                                              separator=input_separator):
                if line_list[scaffold_id_column] not in record_dict_list[-1]:
                    record_dict_list[-1][line_list[scaffold_id_column]] = []
                record_dict_list[-1][line_list[scaffold_id_column]].append([(int(line_list[feature_start_column]) - 1 if coordinates_type == "1-based" else int(line_list[feature_start_column])),
                                                                            int(line_list[feature_end_column])])

        unified_dict = OrderedDict()
        merged_dict = OrderedDict()

        #print record_dict_list[0]

        scaffold_set = set()
        for record_dict in record_dict_list:
            scaffold_set |= set(record_dict.keys())

        for scaffold in scaffold_set:
            unified_dict[scaffold] = []
            merged_dict[scaffold] = []

        for record_dict in record_dict_list:
            for scaffold in record_dict:

                unified_dict[scaffold] += record_dict[scaffold]
                #print "AAAAAAAAAA"
                #print scaffold, unified_dict[scaffold], record_dict[scaffold]

        for scaffold in unified_dict:
            if unified_dict[scaffold]:
                unified_dict[scaffold].sort()
            if unified_dict[scaffold] is None:
                print scaffold

        #print unified_dict

        for scaffold in unified_dict:

            number_of_records = len(unified_dict[scaffold])
            if number_of_records == 0:
                continue

            # [a, b) [c, d), a < b, c < d
            # after sorting c >= a
            i = 1

            prev_coordinates = deepcopy(unified_dict[scaffold][0])

            #print scaffold, number_of_records, prev_coordinates
            #print "\t", unified_dict[scaffold]

            while i < number_of_records:
                if unified_dict[scaffold][i][0] > prev_coordinates[1]: # c > b
                    #print "AAAAAA", "\t", prev_coordinates, unified_dict[scaffold][i]
                    merged_dict[scaffold].append(deepcopy(prev_coordinates))
                    prev_coordinates = deepcopy(unified_dict[scaffold][i])

                elif unified_dict[scaffold][i][1] > prev_coordinates[1]: # d > b; c<=b
                    #print "BBBBBB", "\t",prev_coordinates, unified_dict[scaffold][i]
                    prev_coordinates[1] = deepcopy(unified_dict[scaffold][i][1])
                else: # d <= b
                    #print "CCCCCC", "\t",prev_coordinates, unified_dict[scaffold][i]
                    pass
                i += 1
            if merged_dict[scaffold]:
                if prev_coordinates != merged_dict[scaffold][-1]:
                    merged_dict[scaffold].append(prev_coordinates)
            else:
                merged_dict[scaffold].append(prev_coordinates)

            #print "\t", unified_dict[scaffold]
            #print "\t", merged_dict[scaffold]
        #print unified_dict
        #print merged_dict
        if output_file:
            with self.metaopen(output_file, "w") as out_fd:
                for scaffold in merged_dict:
                    for feature in merged_dict[scaffold]:
                        out_fd.write(output_separator.join(map(str,
                                                               [scaffold,
                                                                feature[0] + 1 if coordinates_type == "1-based" else feature[0],
                                                                feature[1]])) + "\n")

        if return_seqfeature_dict and feature_type:
            feature_dict = OrderedDict()
            for region in merged_dict:
                feature_dict[region] = []
                for (start, stop) in merged_dict[region]:
                    feature_dict[region].append(SeqFeature(FeatureLocation(start, stop),
                                                           type=feature_type,
                                                           strand=None))
            return feature_dict
        elif return_seqfeature_dict and (not feature_type):
            raise ValueError("ERROR!!! Feature type for seqfeature records was not set!")
        else:
            return merged_dict

    def rename_scaffolds_in_gff(self, input_gff, syn_file, output_prefix, verbose=True):

        syn_dict = SynDict(filename=syn_file)
        skipped_id_list = IdSet()

        output_gff = "%s.renamed.gff" % output_prefix
        skipped_gff = "%s.skipped.gff" % output_prefix
        skipped_id_file = "%s.skipped_scaffolds.ids" % output_prefix

        with self.metaopen(input_gff, "r") as in_fd, \
             self.metaopen(output_gff, "w") as out_fd, \
             self.metaopen(skipped_gff, "w") as skipped_fd:

            for line in in_fd:
                if line[0] == "#":
                    out_fd.write(line)
                gff_list = line.split("\t")
                if gff_list[0] in syn_dict:
                    gff_list[0] = syn_dict[gff_list[0]]
                    out_fd.write("\t".join(gff_list))
                else:
                    skipped_fd.write(line)
                    skipped_id_list.add(gff_list[0])

        if verbose:
            print("Not renamed scaffolds: %i" % len(skipped_id_list))

        skipped_id_list.write(skipped_id_file)

    @staticmethod
    def feature_list_entry_to_tab_str(feature_entry):
        return "%s\t%s\t%s\t%s" % (feature_entry[0],
                                   str(feature_entry[1]),
                                   str(feature_entry[2]),
                                   str(feature_entry[3]))

    @staticmethod
    def feature_list_entry_to_gatk_interval_str(feature_entry):
        return "%s:%s-%s" % (feature_entry[0],
                             str(feature_entry[1]),
                             str(feature_entry[2]))

    def get_feature_dict(self, input_gff, output_prefix=None, feature_type_list=["CDS"], unification_key="Parent"):

        feature_dict = SynDict()
        for line_list in self.file_line_as_list_generator(input_gff, comments_prefix="#", separator="\t"):
            annotation_dict = self.parse_gff_annotation_string_to_dict(line_list[self.gff_attribute_column])

            if line_list[self.gff_featuretype_column] not in feature_type_list:
                continue

            if unification_key not in annotation_dict:
                continue

            if annotation_dict[unification_key] not in feature_dict:
                feature_dict[annotation_dict[unification_key]] = []

            feature_dict[annotation_dict[unification_key]].append([line_list[self.gff_scaffold_column],
                                                                   line_list[self.gff_start_column],
                                                                   line_list[self.gff_end_column],
                                                                   line_list[self.gff_strand_column]])

        if output_prefix:
            feature_dict.write("%s.tab" % output_prefix,
                               value_expression=self.feature_list_entry_to_tab_str,
                               line_per_value=True)
            feature_dict.write("%s.coordinates_only.tab" % output_prefix,
                               value_expression=self.feature_list_entry_to_tab_str,
                               line_per_value=True,
                               values_only=True)

            feature_dict.write("%s.list" % output_prefix,
                               value_expression=self.feature_list_entry_to_gatk_interval_str,
                               line_per_value=True)
            feature_dict.write("%s.coordinates_only.list" % output_prefix,
                               value_expression=self.feature_list_entry_to_gatk_interval_str,
                               line_per_value=True,
                               values_only=True)

        return feature_dict







