__author__ = 'mahajrod'
import os
import re
import sys
import math
import pickle

from random import randint

from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

from CustomCollections.GeneralCollections import TwoLvlDict, SynDict, IdList, IdSet
from Routines.File import FileRoutines


class SequenceRoutines(FileRoutines):

    def __init__(self):
        FileRoutines.__init__(self)

        self.ambiguous_nucleotides_dict = OrderedDict({"R": {"A", "G"},
                                                       "Y": {"C", "T"},
                                                       "S": {"C", "G"},
                                                       "W": {"A", "T"},
                                                       "K": {"G", "T"},
                                                       "M": {"A", "C"},
                                                       "B": {"C", "G", "T"},
                                                       "D": {"A", "G", "T"},
                                                       "H": {"A", "C", "T"},
                                                       "V": {"A", "C", "G"},
                                                       "N": {"A", "C", "G", "T"}})

        self.ambiguous_nucleotides_string_dict = OrderedDict({"R": "AG",
                                                              "Y": "CT",
                                                              "S": "CG",
                                                              "W": "AT",
                                                              "K": "GT",
                                                              "M": "AC",
                                                              "B": "CGT",
                                                              "D": "AGT",
                                                              "H": "ACT",
                                                              "V": "ACG",
                                                              "N": "ACGT"})

        self.ambiguous_nucleotides_string_reverse_dict = OrderedDict({"AG": "R",
                                                                      "CT": "Y",
                                                                      "CG": "S",
                                                                      "AT": "W",
                                                                      "GT": "K",
                                                                      "AC": "M",
                                                                      "CGT": "B",
                                                                      "AGT": "D",
                                                                      "ACT": "H",
                                                                      "ACG": "V",
                                                                      "ACGT": "N"})

    def split_fasta(self, input_fasta, output_dir, num_of_recs_per_file=None, num_of_files=None, output_prefix=None,
                    parsing_mode="parse", index_file=None):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        self.safe_mkdir(output_dir)
        out_prefix = self.split_filename(input_fasta)[1] if output_prefix is None else output_prefix

        sequence_dict = self.parse_seq_file(input_fasta, parsing_mode, format="fasta", index_file=index_file) # SeqIO.index_db("temp.idx", input_fasta, "fasta")

        split_index = 0
        records_written = 0
        record_ids_list = list(sequence_dict.keys())
        number_of_records = len(record_ids_list)

        num_of_recs = int(number_of_records/num_of_files) + 1 if num_of_files else num_of_recs_per_file
        while (records_written + num_of_recs) <= number_of_records:
            split_index += 1
            SeqIO.write(self.record_by_id_generator(sequence_dict,
                                                    record_ids_list[records_written:records_written+num_of_recs]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

            records_written += num_of_recs

        if records_written != number_of_records:
            split_index += 1
            SeqIO.write(self.record_by_id_generator(sequence_dict,
                                                    record_ids_list[records_written:]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
        if index_file:
            os.remove(index_file)

        return split_index

    def prepare_region_list_by_length(self, max_length=500000, max_seq_number=10,
                                      length_dict=None, reference=None, parsing_mode="parse", output_dir=None,
                                      split_scaffolds=True, min_scaffold_length=None, black_list_scaffolds=None,
                                      region_file_format='simple'):

        raw_len_dict = length_dict if length_dict else self.get_lengths(record_dict=self.parse_seq_file(reference,
                                                                                                        mode=parsing_mode),
                                                                        out_file=None,
                                                                        close_after_if_file_object=False)
        len_dict = OrderedDict()

        if black_list_scaffolds and min_scaffold_length:
            for scaffold_id in raw_len_dict:
                if (scaffold_id in black_list_scaffolds) or (raw_len_dict[scaffold_id] < min_scaffold_length):
                    continue
                len_dict[scaffold_id] = raw_len_dict[scaffold_id]

        elif black_list_scaffolds:
            for scaffold_id in raw_len_dict:
                if scaffold_id in black_list_scaffolds:
                    continue
                len_dict[scaffold_id] = raw_len_dict[scaffold_id]

        elif min_scaffold_length:
            for scaffold_id in raw_len_dict:
                if raw_len_dict[scaffold_id] < min_scaffold_length:
                    continue
                len_dict[scaffold_id] = raw_len_dict[scaffold_id]

        else:
            len_dict = raw_len_dict

        scaffold_ids = IdList(len_dict)
        scaffold_ids.write("%s/scaffold.ids" % output_dir)
        len_dict.write("%s/scaffold.len" % output_dir)

        number_of_scaffolds = len(len_dict)
        max_length_soft_threshold = None if max_length is None else int(1.5 * max_length)
        region_list = []

        remnant_seq_list = []
        remnant_seq_length = 0

        scaffold_to_region_correspondence_dict = SynDict()
        region_index = 0

        if max_length_soft_threshold is None:
            key_list = list(len_dict.keys())
            bins = np.arange(0, number_of_scaffolds, max_seq_number)
            bins = bins if bins[-1] == number_of_scaffolds else np.append(bins, number_of_scaffolds)

            for i in range(0, len(bins)-1):
                region_list.append(key_list[bins[i]:bins[i+1]])
                for scaffold in key_list[bins[i]:bins[i+1]]:
                    scaffold_to_region_correspondence_dict[scaffold] = [i]

        elif not split_scaffolds:
            bunch_length = 0
            bunch_list = []

            for region in len_dict:
                if len_dict[region] >= max_length:
                    region_list.append([region])
                    scaffold_to_region_correspondence_dict[region] = [region_index]
                    region_index += 1
                else:
                    bunch_list.append(region)
                    bunch_length += len_dict[region]
                    if (bunch_length >= max_length) or (len(bunch_list) == max_seq_number):
                        region_list.append(bunch_list)
                        for scaffold in bunch_list:
                            scaffold_to_region_correspondence_dict[scaffold] = [region_index]
                        region_index += 1
                        bunch_length = 0
                        bunch_list = []
            if bunch_list:
                region_list.append(bunch_list)
                for scaffold in bunch_list:
                    scaffold_to_region_correspondence_dict[scaffold] = [region_index]
                region_index += 1
                bunch_length = 0
                bunch_list = []
        else:
            remnant_seq_length = 0
            for region in len_dict:
                print region
                if len(remnant_seq_list) == max_seq_number:
                    #print region_list
                    #print remnant_seq_list
                    region_list.append(remnant_seq_list)
                    for remnant_entry in remnant_seq_list:
                        if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                        else:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                    region_index += 1
                    remnant_seq_list = []
                    remnant_seq_length = 0

                if len_dict[region] > max_length:
                    points = np.arange(0, len_dict[region], max_length)
                    if len(points) > 1:
                        scaffold_to_region_correspondence_dict[region] = list(np.arange(region_index, region_index + len(points) - 1))
                        region_index += len(points) - 1
                        for i in range(0, len(points) - 1):
                            region_list.append([[region, points[i] + 1, points[i+1]]])

                        #remnant = [region, points[-1] + 1, len_dict[region]]
                        remnant_length = len_dict[region] - points[-1]
                        if remnant_length + max_length <= max_length_soft_threshold:
                            region_list[-1][0][2] = len_dict[region]
                            remnant_length = 0
                        else:
                            region_list.append([[region, points[-1] + 1, len_dict[region]]])
                            scaffold_to_region_correspondence_dict[region].append(region_index)
                            region_index += 1
                    else:
                        #remnant = [region, 1, len_dict[region]]
                        #remnant_length = len_dict[region]
                        region_list.append([[region, 1, len_dict[region]]])
                        scaffold_to_region_correspondence_dict[region] = [region_index]
                        region_index += 1
                    continue
                else:
                    print "BBBBB"

                    remnant = [region, 1, len_dict[region]]
                    remnant_length = len_dict[region]
                    print remnant, remnant_length
                    print "\n"

                if remnant_seq_length + remnant_length <= max_length_soft_threshold:
                    print "AAAAA"
                    print remnant_seq_length, remnant_length, max_length_soft_threshold
                    remnant_seq_list.append(remnant)
                    remnant_seq_length += remnant_length
                    print remnant_seq_list
                    print "\n"
                else:
                    region_list.append(remnant_seq_list)
                    for remnant_entry in remnant_seq_list:
                        #print scaffold
                        if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                        else:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                    region_index += 1
                    print "CCCCCCCCC"
                    remnant_seq_list = [remnant]
                    remnant_seq_length = remnant_length
                    print remnant_seq_list
                    print remnant_seq_length
                    print "\n"
            else:
                if remnant_seq_list:
                    region_list.append(remnant_seq_list)
                    for remnant_entry in remnant_seq_list:
                        if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                        else:
                            scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                    region_index += 1

        if output_dir:
            for directory in output_dir, "%s/splited/" % output_dir:
                self.safe_mkdir(directory)

            index = 1
            for regions in region_list:
                with open("%s/splited/region_%i.list" % (output_dir, index), "w") as out_fd:
                    for region in regions:
                        if isinstance(region, str):
                            out_fd.write(region)
                            out_fd.write("\n")
                        else:
                            if region_file_format == 'simple':
                                if len(region) == 3:
                                    out_fd.write("%s\t%s\t%s\n" % (region[0], region[1], region[2]))
                                elif len(region) == 1:
                                    out_fd.write(region[0])
                                    out_fd.write("\n")
                            elif region_file_format == 'GATK':
                                if len(region) == 3:
                                    out_fd.write("%s:%s-%s\n" % (region[0], region[1], region[2]))
                                elif len(region) == 1:
                                    out_fd.write(region[0])
                                    out_fd.write("\n")

                index += 1
            scaffold_to_region_correspondence_dict.write("%s/SCAFFOLD_TO_REGION.correspondence" % output_dir,
                                                         splited_values=True)
        return region_list, scaffold_to_region_correspondence_dict

    def split_fasta_by_seq_len(self, input_fasta, output_dir, max_len_per_file=None, output_prefix=None,
                               parsing_mode="parse", index_file='temp.idx'):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        self.safe_mkdir(output_dir)

        out_prefix = self.split_filename(input_fasta)[1] if output_prefix is None else output_prefix
        sequence_dict = self.parse_seq_file(input_fasta, parsing_mode, format="fasta", index_file=index_file)
        length = 0

        for record_id in sequence_dict:
            length += len(sequence_dict[record_id].seq)

        max_len = max_len_per_file if max_len_per_file else int(length / self.threads)

        split_index = 1
        id_list = []
        total_length = 0

        for record_id in sequence_dict:
            record_length = len(sequence_dict[record_id].seq)
            if record_length >= max_len:
                SeqIO.write(sequence_dict[record_id],
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

            elif total_length + record_length > max_len:
                SeqIO.write(self.record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = record_length
                id_list = [record_id]

            elif total_length + record_length == max_len:
                id_list.append(record_id)
                SeqIO.write(self.record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = 0
                id_list = []

            elif total_length + record_length < max_len:
                id_list.append(record_id)
                total_length += record_length
                continue

            split_index += 1

        if id_list:
            SeqIO.write(self.record_by_id_generator(sequence_dict, id_list),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
        if parsing_mode == "index_db":
            os.remove(index_file)

    def extract_common_sequences(self, list_of_files_with_sequences_of_samples, list_of_names_of_samples,
                                 output_dir, separator="_", format="fasta"):

        def generator_with_id_correction(samples_seq_dict, common_record_id):
            for sample in samples_seq_dict:
                record = samples_seq_dict[sample][common_record_id]
                record.id = "%s%s%s" % (sample, separator, record.id)
                yield record

        self.safe_mkdir(output_dir)
        index = 0
        samples_seq_dict = OrderedDict()
        for filename, sample_name in zip(list_of_files_with_sequences_of_samples, list_of_names_of_samples):
            samples_seq_dict[sample_name] = SeqIO.index_db("tmp_%i.idx" % index, filename, format=format)
            index += 1
        common_sequence_ids = set(samples_seq_dict[list_of_names_of_samples[0]].keys())

        for sample_name in list_of_names_of_samples[1:]:
            common_sequence_ids = common_sequence_ids & set(samples_seq_dict[sample_name].keys())

        for common_id in common_sequence_ids:
            SeqIO.write(generator_with_id_correction(samples_seq_dict, common_id),
                        "%s%s.%s" % (self.check_path(output_dir), common_id, format),
                        format=format)
        for i in range(0, index):
            os.remove("tmp_%i.idx" % i)

    @staticmethod
    def get_cds_by_coordinates_from_transcript_generator(transcript_dict, coordinates_dict, zero_based=False):
        for record_id in transcript_dict:
            if record_id not in coordinates_dict:
                print("WARNING!!! No CDS coordinates for transcript %s" % record_id)
                continue

            start = coordinates_dict[record_id][0] if zero_based else coordinates_dict[record_id][0] - 1
            end = coordinates_dict[record_id][1] + 1 if zero_based else coordinates_dict[record_id][1]
            new_record = deepcopy(transcript_dict[record_id])
            new_record.seq = new_record.seq[start:end]
            yield new_record

    def get_cds_by_bed_from_transcript_file(self, transcript_file, coordinate_bed, output_file, zero_based=False,
                                            transcript_format="fasta", parsing_mode="parse"):
        transcript_dict = self.parse_seq_file(transcript_file, format=transcript_format,
                                              mode=parsing_mode)
        coordinates_dict = OrderedDict()
        with open(coordinate_bed, 'r') as in_fd:
            for line in in_fd:
                line_list = line.strip().split("\t")
                coordinates_dict[line_list[0]] = (int(line_list[1]), int(line_list[2]))

        SeqIO.write(self.get_cds_by_coordinates_from_transcript_generator(transcript_dict,
                                                                          coordinates_dict,
                                                                          zero_based=zero_based),
                    output_file,
                    transcript_format)

    @staticmethod
    def get_lengths(record_dict, out_file=None, close_after_if_file_object=False):
        #engths_dict = SynDict()
        length_list = []
        for record_id in record_dict:
            length_list.append((record_id, len(record_dict[record_id])))
            #lengths_dict[record_id] = len(record_dict[record_id])

        length_list.sort(key=lambda s: s[1], reverse=True)
        lengths_dict = SynDict(length_list)
        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)

        return lengths_dict

    @staticmethod
    def get_first_seq_length_from_file(seq_file, format="fasta"):
        for record in SeqIO.parse(seq_file, format=format):
            return len(record.seq)

    @staticmethod
    def get_first_record_from_file(seq_file, format="fasta"):
        for record in SeqIO.parse(seq_file, format=format):
            return record.seq

    def get_lengths_from_seq_file(self, input_file_list, format="fasta", out_file=None, close_after_if_file_object=False,
                                  parsing_mode="parse"):
        record_dict = self.parse_seq_file(input_file_list, mode=parsing_mode, format=format, index_file="tmp.idx")

        lengths_dict = SynDict()
        for record_id in record_dict:
            lengths_dict[record_id] = len(record_dict[record_id])

        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)
        if parsing_mode == "index_db":
            os.remove("tmp.idx")
        return lengths_dict

    @staticmethod
    def record_by_id_generator(record_dict, id_list=[], verbose=False, coincidence_mode="exact",
                               allow_multiple_coincidence_report=False, invert_match=False):
        if id_list:
            if invert_match:
                for record_id in record_dict:
                    if record_id not in id_list:
                        yield record_dict[record_id]
            else:
                for record_id in id_list:
                    if coincidence_mode == "exact":
                        if record_id in record_dict:
                            yield record_dict[record_id]
                        else:
                            if verbose:
                                sys.stderr.write("Not found: %s\n" % record_id)
                    elif coincidence_mode == "partial":
                        #print "AAAAAA"
                        entry_list = []
                        if record_id in record_dict:
                            #print "BBBBBB"
                            yield record_dict[record_id]
                        else:
                            #print "CCCCCCCCCC"
                            for dict_entry in record_dict:
                                if record_id in dict_entry:
                                    entry_list.append(dict_entry)
                            if len(entry_list) > 1:
                                if allow_multiple_coincidence_report:
                                    sys.stderr.write("WARNING!!! Multiple coincidence for %s" % record_id)
                                    sys.stderr.write("\treporting all...")
                                else:
                                    sys.stderr.write("ERROR!!! Multiple coincidence for %s" % record_id)
                                    raise ValueError("Multiple coincidence for %s" % record_id)
                               # print entry_list
                                for entry in entry_list:
                                    yield record_dict[entry]
                            elif len(entry_list) == 1:
                                yield record_dict[entry_list[0]]
                        if (not entry_list) and verbose:
                            sys.stderr.write("Not found: %s\n" % record_id)
                    else:
                        raise ValueError("Unknown coincidence mode: %s" % coincidence_mode)
        else:
            for record_id in record_dict:
                yield record_dict[record_id]

    @staticmethod
    def record_from_dict_generator(record_dict):
        for record_id in record_dict:
            yield record_dict[record_id]

    def extract_sequence_by_ids(self, sequence_file, id_file, output_file, format="fasta", verbose=False,
                                id_column_number=0, coincidence_mode="exact", allow_multiple_coincidence_report=False,
                                syn_file=None, parsing_mode="parse", index_file="tmp.idx", invert_match=False):

        id_list = IdList()
        id_list.read(id_file, column_number=id_column_number)

        converted_ids = IdList()
        if syn_file:
            syn_dict = SynDict(filename=syn_file)
            for entry in id_list:
                if entry not in syn_dict:
                    if verbose:
                        print("No synonym for %s" % entry)
                    converted_ids.append(entry)
                    continue
                converted_ids.append(syn_dict[entry])
            id_list = converted_ids

        if verbose:
            print("Parsing %s..." % (sequence_file if isinstance(id_file, str) else ",".join(id_file)))

        sequence_dict = self.parse_seq_file(sequence_file,parsing_mode, format=format, index_file=index_file) # SeqIO.index_db(tmp_index_file, sequence_file, format=format)
        SeqIO.write(self.record_by_id_generator(sequence_dict, id_list, verbose=verbose,
                                                coincidence_mode=coincidence_mode,
                                                allow_multiple_coincidence_report=allow_multiple_coincidence_report,
                                                invert_match=invert_match),
                    output_file, format=format)
        if parsing_mode == "index_db":
            os.remove(index_file)

    def extract_sequences_by_length_from_file(self, input_file, output_file, min_len=1, max_len=None, format="fasta",
                                              tmp_index_file="tmp.idx", id_file=None, parsing_mode='parse'):

        if (min_len is None) and (max_len is None):
            raise ValueError("Both minimum and maximum lengths were not set")
        elif (min_len is not None) and (max_len is not None) and (min_len > max_len):
            raise ValueError("Minimum length is greater then maximum lengths")
        elif (min_len is not None) and (min_len < 0):
            raise ValueError("Minimum length is below zero")
        elif (max_len is not None) and (max_len < 0):
            raise ValueError("Maximum length is below zero")

        sequence_dict = self.parse_seq_file(input_file, mode=parsing_mode, format=format, index_file=tmp_index_file) #        SeqIO.index_db(tmp_index_file, input_file, format=format)

        if (min_len is not None) and (min_len > 1) and (max_len is not None):
            length_expression = lambda record: min_len <= len(record.seq) <= max_len
        elif (min_len is not None) and (min_len > 1):
            length_expression = lambda record: len(record.seq) >= min_len
        elif max_len is not None:
            length_expression = lambda record: len(record.seq) <= max_len
        else:
            length_expression = lambda record: True

        SeqIO.write(self.record_by_expression_generator(sequence_dict, expression=length_expression,
                                                        id_file=id_file),
                    output_file, format=format)
        if parsing_mode == "index_db":
            os.remove(tmp_index_file)

    @staticmethod
    def find_gaps(record_dict):
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        gaps_dict = {}
        for region in record_dict:
            gaps_dict[region] = SeqRecord(seq=record_dict[region].seq,
                                          id=record_dict[region].id,
                                          description=record_dict[region].description)
            gaps = gap_reg_exp.finditer(str(record_dict[region].seq))  # iterator with
            for match in gaps:
                gaps_dict[region].features.append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                  type="gap", strand=None))
        return gaps_dict

    @staticmethod
    def record_by_expression_generator(record_dict, expression=None, id_file="passed_records.ids"):
        """
        :param record_dict: dictionary containing Biopython SeqRecords as values
        :param expression: function to apply to all records in record_dict. If it returns True record will be yielded
        :param id_file: file to write ids of records passed expression
        :return: None
        """
        if id_file:
            with open(id_file, "w") as id_fd:
                if expression is None:
                    for record_id in record_dict:
                        id_fd.write(record_id + "\n")
                        yield record_dict[record_id]
                else:
                    for record_id in record_dict:
                        if expression(record_dict[record_id]):
                            id_fd.write(record_id + "\n")
                            yield record_dict[record_id]
        else:
            if expression is None:
                for record_id in record_dict:
                    yield record_dict[record_id]
            else:
                for record_id in record_dict:
                    if expression(record_dict[record_id]):
                        yield record_dict[record_id]

    @staticmethod
    def filter_seq_by_expression(record_dict, expression):
        filtered_ids = IdSet()
        filtered_out_ids = IdSet()
        for record_id in record_dict:
            if expression(record_dict[record_id]):
                filtered_ids.add(record_id)
            else:
                filtered_out_ids.add(record_id)
        return filtered_ids, filtered_out_ids

    @staticmethod
    def expression_for_seq_ids(record, compiled_reg_expression):
        return True if compiled_reg_expression.search(record.id)else False

    def filter_seq_by_ids_and_compiled_reg_expression(self, record_dict, compiled_reg_expression):

        def expression(record):
            return True if compiled_reg_expression.search(record.id)else False

        return self.filter_seq_by_expression(record_dict, expression)

    def filter_seq_by_ids_and_reg_expression(self, record_dict, reg_expression, reg_exp_flags=0):

        compiled_reg_expression = re.compile(reg_expression, flags=reg_exp_flags)

        def expression(record):
            return True if compiled_reg_expression.search(record.id)else False

        return self.filter_seq_by_expression(record_dict, expression)

    def filter_seq_by_expression_from_file(self, input_file, expression_function,
                                           output_filtered_file, output_filtered_out_file,
                                           parsing_mode="index_db", format="fasta",
                                           index_file="tmp.idx", retain_index=False):
        record_dict = self.parse_seq_file(input_file, parsing_mode, format=format, index_file=index_file)

        filtered_ids, filtered_out_ids = self.filter_seq_by_expression(record_dict, expression_function)

        SeqIO.write(self.record_by_id_generator(record_dict, filtered_ids), output_filtered_file, format=format)
        SeqIO.write(self.record_by_id_generator(record_dict, filtered_out_ids), output_filtered_out_file, format=format)

        if (parsing_mode == "index_db") and (not retain_index):
            os.remove(index_file)

    def filter_seq_by_reg_expression_from_file(self, input_file, reg_expression,
                                               output_filtered_file, output_filtered_out_file,
                                               parsing_mode="index_db", format="fasta",
                                               index_file="tmp.idx", retain_index=False, reg_exp_flags=0):

        compiled_reg_expression = re.compile(reg_expression, flags=reg_exp_flags)

        def expression(record):
            return True if compiled_reg_expression.search(record.id)else False

        self.filter_seq_by_expression_from_file(input_file, expression, output_filtered_file, output_filtered_out_file,
                                                parsing_mode=parsing_mode, format=format,
                                                index_file=index_file, retain_index=retain_index)

    def check_pairing(self, forward_record_dict, reverse_record_dict, forward_record_id_suffix,
                      reverse_record_id_suffix):
        forward_suffix_len = len(forward_record_id_suffix)
        reverse_suffix_len = len(reverse_record_id_suffix)

        forward_ids = forward_record_dict.keys()
        forward_ids_common_part = []
        for forward_id in forward_ids:
            if forward_id[-forward_suffix_len:] != forward_record_id_suffix:
                raise ValueError("Forward record %s doesn't have forward id suffix" % forward_id)
            forward_ids_common_part.append(forward_id[:-forward_suffix_len])

        reverse_ids = reverse_record_dict.keys()
        reverse_ids_common_part = []
        for reverse_id in reverse_ids:
            if reverse_id[-reverse_suffix_len:] != reverse_record_id_suffix:
                raise ValueError("Reverse record %s doesn't have forward id suffix" % reverse_id)
            reverse_ids_common_part.append(reverse_id[:-reverse_suffix_len])

        paired_ids = self.intersect_ids([forward_ids_common_part], [reverse_ids_common_part], mode="common")
        only_forward_ids = self.intersect_ids(forward_ids_common_part, reverse_ids_common_part, mode="only_a")
        only_reverse_ids = self.intersect_ids(forward_ids_common_part, reverse_ids_common_part, mode="only_b")

        forward_full_paired_ids = [common_id + forward_record_id_suffix for common_id in paired_ids]
        reverse_full_paired_ids = [common_id + reverse_record_id_suffix for common_id in paired_ids]
        only_forward_full_ids = [only_forward_id + forward_record_id_suffix for only_forward_id in only_forward_ids]
        only_reverse_full_ids = [only_reverse_id + reverse_record_id_suffix for only_reverse_id in only_reverse_ids]

        return forward_full_paired_ids, reverse_full_paired_ids, only_forward_full_ids, only_reverse_full_ids

    def check_pairing_from_file(self, forward_file, reverse_file, output_prefix, forward_record_id_suffix,
                                reverse_record_id_suffix, parsing_mode="index_db", format="fasta",
                                forward_index_file="forward_tmp.idx", reverse_index_file="reverse_tmp.idx",
                                retain_index=False, output_file_extension="fasta"):
        forward_dict = self.parse_seq_file(forward_file, parsing_mode, format=format, index_file=forward_index_file)
        reverse_dict = self.parse_seq_file(reverse_file, parsing_mode, format=format, index_file=reverse_index_file)

        forward_full_paired_ids, \
        reverse_full_paired_ids, \
        only_forward_full_ids, \
        only_reverse_full_ids = self.check_pairing(forward_dict, reverse_dict,
                                                   forward_record_id_suffix, reverse_record_id_suffix)

        forward_paired_file = "%s.pe.forward.%s" % (output_prefix, output_file_extension)
        reverse_paired_file = "%s.pe.reverse.%s" % (output_prefix, output_file_extension)
        forward_unpaired_file = "%s.se.forward.%s" % (output_prefix, output_file_extension)
        reverse_unpaired_file = "%s.se.reverse.%s" % (output_prefix, output_file_extension)

        for (dictionary, ids, filename) in zip((forward_dict, reverse_dict, forward_dict, reverse_dict),
                                               (forward_full_paired_ids, reverse_full_paired_ids, only_forward_full_ids, only_reverse_full_ids),
                                               (forward_paired_file, reverse_paired_file, forward_unpaired_file, reverse_unpaired_file)):
            SeqIO.write(self.record_by_id_generator(dictionary, ids, verbose=True), filename, format=format)

        if (parsing_mode == "index_db") and (not retain_index):
            os.remove(forward_index_file)
            os.remove(reverse_index_file)
    """
    @staticmethod
    def check_pairing_from_interleaved_dict(record_dict, forward_record_id_suffix, reverse_record_id_suffix):
        forward_suffix_len = len(forward_record_id_suffix)
        reverse_suffix_len = len(reverse_record_id_suffix)
    """

    @staticmethod
    def parse_seq_file(input_file, mode, format="fasta", index_file=None):
        if mode == "index_db" or ((not isinstance(input_file, str)) and (len(input_file) > 1)):
            index = index_file if index_file else "tmp.idx"
            seq_dict = SeqIO.index_db(index, [input_file] if isinstance(input_file, str) else input_file, format=format)
        elif mode == "index":
            seq_dict = SeqIO.index(input_file if isinstance(input_file, str) else input_file[0], format=format)
        elif mode == "parse":
            seq_dict = OrderedDict()
            for record in SeqIO.parse(input_file if isinstance(input_file, str) else input_file[0], format=format):
                seq_dict[record.id] = record
            #seq_dict = SeqIO.to_dict(SeqIO.parse(input_file if isinstance(input_file, str) else input_file[0], format=format))

        return seq_dict

    def translated_seq_generator(self, cds_dict, id_expression=None, genetic_code_table=1, translate_to_stop=True,
                                 stop_codon_symbol="*", prefix_of_file_inframe_stop_codons_seqs=None):

        seqs_inframe_stop_codons_ids = IdSet()
        for cds_id in cds_dict:
            pep_seq = cds_dict[cds_id].seq.translate(to_stop=translate_to_stop, table=genetic_code_table,
                                                     stop_symbol=stop_codon_symbol)
            pep_id = id_expression(cds_id) if id_expression else cds_id
            description = " cds_id=%s" % cds_id
            pep_record = SeqRecord(id=pep_id, description=description, seq=pep_seq)

            if stop_codon_symbol in pep_seq:
                seqs_inframe_stop_codons_ids.add(pep_id)

            yield pep_record

        if prefix_of_file_inframe_stop_codons_seqs:
            id_file = "%s.ids" % prefix_of_file_inframe_stop_codons_seqs
            seq_file = "%s.cds" % prefix_of_file_inframe_stop_codons_seqs

            seqs_inframe_stop_codons_ids.write(id_file)
            SeqIO.write(self.record_by_id_generator(cds_dict, seqs_inframe_stop_codons_ids), seq_file,
                        format="fasta")

    def translate_sequences_from_file(self, input_file, output_file, format="fasta", id_expression=None,
                                      genetic_code_table=1, translate_to_stop=True,
                                      prefix_of_file_inframe_stop_codons_seqs="cds_with_in_frame_stop_codons",
                                      mode="parse", index_file="tmp.idx"):

        cds_dict = self.parse_seq_file(input_file, mode, format=format, index_file=index_file)
        #cds_dict = SeqIO.index_db("tmp.idx", input_file, format=format)
        SeqIO.write(self.translated_seq_generator(cds_dict, id_expression=id_expression,
                                                  genetic_code_table=genetic_code_table,
                                                  translate_to_stop=translate_to_stop,
                                                  prefix_of_file_inframe_stop_codons_seqs=prefix_of_file_inframe_stop_codons_seqs),
                    output_file, format=format)
        if mode == "index_db":
            os.remove("tmp.idx")

    @staticmethod
    def compare_sequences(record_dict_1, record_dict_2):
        id_set_dict_1 = IdSet(record_dict_1.keys())
        id_set_dict_2 = IdSet(record_dict_2.keys())
        common_ids = id_set_dict_1 & id_set_dict_2
        id_set_unique_to_dict_1 = id_set_dict_1 - id_set_dict_2
        id_set_unique_to_dict_2 = id_set_dict_2 - id_set_dict_1

        equal_seq_ids_set = IdSet()
        non_equal_seq_ids_set = IdSet()
        for seq_id in common_ids:
            if record_dict_1[seq_id].seq == record_dict_2[seq_id].seq:
                equal_seq_ids_set.add(seq_id)
            else:
                non_equal_seq_ids_set.add(seq_id)

        return equal_seq_ids_set, non_equal_seq_ids_set, id_set_unique_to_dict_1, id_set_unique_to_dict_2

    def compare_sequences_from_files(self, seq_file_1, seq_file_2, output_prefix, format="fasta", verbose=True):
        record_dict_1 = SeqIO.index_db("tmp_A.idx", seq_file_1, format=format)
        record_dict_2 = SeqIO.index_db("tmp_B.idx", seq_file_2, format=format)

        comparison_results = self.compare_sequences(record_dict_1, record_dict_2)
        number_of_same_seqs = len(comparison_results[0])
        number_of_different_seqs = len(comparison_results[1])
        number_of_unique_to_file_1_seqs = len(comparison_results[2])
        number_of_unique_to_file_2_seqs = len(comparison_results[3])
        for i, extension in zip(range(0, 4), ("same", "different", "unique_to_file_A", "unique_to_file_B")):
            comparison_results[i].write("%s.%s" % (output_prefix, extension))

        for filename in "tmp_A.idx", "tmp_B.idx":
            os.remove(filename)
        stat_string = "Same sequences: %i\nDifferent sequences: %i\nUnique to file A: %i\nUnique to file B : %i\n" % \
                      (number_of_same_seqs, number_of_different_seqs,
                       number_of_unique_to_file_1_seqs, number_of_unique_to_file_2_seqs)

        with open("%s.stat" % output_prefix, "w") as stat_fd:
            stat_fd.write(stat_string)

        if verbose:
            print(stat_string)

        return comparison_results

    @staticmethod
    def check_proteins_for_stop_codons(protein_dict, stop_codon_symbol_set=("*", ".")):
        pep_with_stop_codons_ids = IdSet()
        for pep_id in protein_dict:
            for stop_codon_symbol in stop_codon_symbol_set:
                if stop_codon_symbol in protein_dict[pep_id].seq:
                    pep_with_stop_codons_ids.add(pep_id)
                    continue
        return pep_with_stop_codons_ids

    def check_proteins_for_stop_codons_from_file(self, pep_file, output_prefix, stop_codon_symbol_set=("*", "."), format="fasta",
                                                 parsing_mode="parse"):
        pep_with_stop_codons_ids_file = "%s.ids" % output_prefix
        pep_with_stop_codons_seq_file = "%s.pep" % output_prefix
        pep_dict =  self.parse_seq_file(pep_file, parsing_mode, format=format, index_file="tmp.idx")

        pep_with_stop_codons_ids = self.check_proteins_for_stop_codons(pep_dict,
                                                                       stop_codon_symbol_set=stop_codon_symbol_set)
        pep_with_stop_codons_ids.write(pep_with_stop_codons_ids_file)
        SeqIO.write(self.record_by_id_generator(pep_dict, pep_with_stop_codons_ids), pep_with_stop_codons_seq_file,
                    format=format)
        if parsing_mode == "index_db":
            os.remove("tmp.idx")
        return pep_with_stop_codons_ids

    @staticmethod
    def trim_proteins_by_stop(record_dict, stop_codon_symbol_set=("*", ".")):
        trimmed_record_dict = OrderedDict()
        trimmed_ids = IdList()
        reg_exp = "|".join(map(lambda symbol: "\\" + symbol, stop_codon_symbol_set))
        #print reg_exp
        stop_reg_exp = re.compile(reg_exp)
        for record_id in record_dict:
            stop_match = stop_reg_exp.search(str(record_dict[record_id].seq))

            if stop_match:
                #print stop_match.start
                trimmed_record_dict[record_id] = SeqRecord(id=record_dict[record_id].id,
                                                           description=record_dict[record_id].description,
                                                           seq=record_dict[record_id].seq[:stop_match.start()])
                trimmed_ids.append(record_id)
            else:
                trimmed_record_dict[record_id] = record_dict[record_id]
        return trimmed_record_dict, trimmed_ids

    def trim_proteins_by_stop_from_file(self, pep_file, output_prefix, stop_codon_symbol_set=("*", "."), format="fasta",
                                        parsing_mode="parse"):
        trimmed_pep_id_file = "%s.trimmed.ids" % output_prefix
        trimmed_pep_file = "%s.pep" % output_prefix
        pep_dict = self.parse_seq_file(pep_file, parsing_mode, format=format, index_file="tmp.idx")

        trimmed_pep_dict, trimmed_pep_id_list = self. trim_proteins_by_stop(pep_dict, stop_codon_symbol_set)
        trimmed_pep_id_list.write(trimmed_pep_id_file)
        SeqIO.write(self.record_by_id_generator(trimmed_pep_dict), trimmed_pep_file,
                    format=format)
        if parsing_mode == "index_db":
            os.remove("tmp.idx")
        return trimmed_pep_id_list

    @staticmethod
    def check_cds_for_in_frame_stop_codons(cds_dict, genetic_code_table=1):
        cds_with_stop_codons_ids = IdSet()
        stop_codon_symbol = "*"
        for cds_id in cds_dict:
            pep_seq = cds_dict[cds_id].seq.translate(to_stop=False, table=genetic_code_table,
                                                     stop_symbol=stop_codon_symbol)
            if stop_codon_symbol in pep_seq:
                cds_with_stop_codons_ids.add(cds_id)
                continue
        return cds_with_stop_codons_ids

    def check_cds_for_stop_codons_from_file(self, cds_file, output_prefix, genetic_code_table=1, format="fasta"):
        cds_with_stop_codons_ids_file = "%s.ids" % output_prefix
        cds_with_stop_codons_seq_file = "%s.cds" % output_prefix
        cds_dict = SeqIO.index_db("tmp.idx", cds_file, format=format)

        cds_with_stop_codons_ids = self.check_cds_for_in_frame_stop_codons(cds_dict,
                                                                           genetic_code_table=genetic_code_table)
        cds_with_stop_codons_ids.write(cds_with_stop_codons_ids_file)
        SeqIO.write(self.record_by_id_generator(cds_dict, cds_with_stop_codons_ids), cds_with_stop_codons_seq_file, 
                    format=format)
        os.remove("tmp.idx")
        return cds_with_stop_codons_ids

    @staticmethod
    def check_for_selenocystein_presence(pep_dict):
        selenocystein_ids = IdSet()

        for pep_id in pep_dict:
            if "U" in pep_dict[pep_id].seq:
                selenocystein_ids.add(pep_id)

        return selenocystein_ids

    def check_selenocystein_presence_from_file(self, pep_file, output_prefix, format="fasta"):
        selenocystein_ids_file = "%s.ids" % output_prefix
        selenocystein_pep_file = "%s.pep" % output_prefix
        pep_dict = SeqIO.index_db("tmp.idx", pep_file, format=format)

        selenocystein_ids = self.check_for_selenocystein_presence(pep_dict)
        selenocystein_ids.write(selenocystein_ids_file)
        SeqIO.write(self.record_by_id_generator(pep_dict, selenocystein_ids), selenocystein_pep_file, format=format)
        os.remove("tmp.idx")
        return selenocystein_ids

    @staticmethod
    def get_cds_to_pep_accordance(cds_dict, pep_dict, verbose=False,
                                  parsing_mode="index_db", genetic_code_table=1,
                                  include_id_check=False):
        cds_pep_accordance_dict = SynDict()
        cds_with_absent_proteins_id_list = IdList()

        for cds_id in cds_dict:
            cds_pep = cds_dict[cds_id].seq.translate(to_stop=True, table=genetic_code_table)
            for pep_id in pep_dict:
                if include_id_check:
                    if cds_id in pep_dict:
                        cds_pep_accordance_dict[cds_id] = cds_id
                        if parsing_mode == "parse":
                            pep_dict.pop(pep_id, None)
                            break
                if len(pep_dict[pep_id].seq)*3 > len(cds_dict[cds_id]):
                    continue

                if cds_pep == pep_dict[pep_id].seq:
                    cds_pep_accordance_dict[cds_id] = pep_id
                    if parsing_mode == "parse":
                        pep_dict.pop(pep_id, None)
                    break
            else:
                cds_with_absent_proteins_id_list.append(cds_id)
                if verbose:
                    print("Protein was not found for %s CDS" % cds_id)

        return cds_pep_accordance_dict, cds_with_absent_proteins_id_list

    @staticmethod
    def get_transcript_to_pep_accordance(transcript_dict, pep_dict, verbose=False,
                                         parsing_mode="parse", genetic_code_table=1,
                                         include_id_check=False):

        transcript_pep_accordance_dict = SynDict()
        transcript_with_absent_proteins_id_list = IdList()
        transcript_with_several_proteins_id_list = IdList()

        for transcript_id in transcript_dict:
            transcript_pep_0 = transcript_dict[transcript_id].seq.translate(table=genetic_code_table)
            transcript_pep_1 = transcript_dict[transcript_id].seq[1:].translate(table=genetic_code_table)
            transcript_pep_2 = transcript_dict[transcript_id].seq[2:].translate(table=genetic_code_table)
            for pep_id in pep_dict:
                if include_id_check:
                    if transcript_id in pep_dict:
                        if transcript_id not in transcript_pep_accordance_dict:
                            transcript_pep_accordance_dict[transcript_id] = [transcript_id]
                        else:
                            transcript_pep_accordance_dict[transcript_id].append(transcript_id)
                        if parsing_mode == "parse":
                            pep_dict.pop(pep_id, None)
                            break
                if len(pep_dict[pep_id].seq)*3 > len(transcript_dict[transcript_id]):
                    continue
                if (pep_dict[pep_id].seq in transcript_pep_0) or (pep_dict[pep_id].seq in transcript_pep_1) or (pep_dict[pep_id].seq in transcript_pep_2):
                    if transcript_id not in transcript_pep_accordance_dict:
                        transcript_pep_accordance_dict[transcript_id] = [pep_id]
                    else:
                        transcript_pep_accordance_dict[transcript_id].append(pep_id)

                    if parsing_mode == "parse":
                        pep_dict.pop(pep_id, None)
                    break
            else:
                transcript_with_absent_proteins_id_list.append(transcript_id)
                if verbose:
                    print("Protein was not found for %s transcript" % transcript_id)

        for transcript_id in transcript_pep_accordance_dict:
            if len(transcript_pep_accordance_dict[transcript_id]) > 1:
                transcript_with_several_proteins_id_list.append(transcript_id)

        return transcript_pep_accordance_dict, transcript_with_absent_proteins_id_list, transcript_with_several_proteins_id_list

    def get_cds_to_pep_accordance_from_files(self, cds_file, pep_file, output_file, cds_with_no_pep_idfile=None, format="fasta",
                                             verbose=True, parsing_mode="parse", index_file_suffix="tmp.idx",
                                             genetic_code_table=1, include_id_check=False):

        cds_dict = self.parse_seq_file(cds_file, mode=parsing_mode, format=format,
                                       index_file="transcript_%s" % index_file_suffix)
        pep_dict = self.parse_seq_file(pep_file, mode=parsing_mode, format=format,
                                       index_file="pep_%s" % index_file_suffix)

        cds_pep_accordance_dict, cds_with_absent_proteins_id_list = self.get_cds_to_pep_accordance(cds_dict, pep_dict, verbose=verbose,
                                                                                                   parsing_mode=parsing_mode,
                                                                                                   genetic_code_table=genetic_code_table,
                                                                                                   include_id_check=include_id_check)
        cds_pep_accordance_dict.write(output_file)
        if cds_with_no_pep_idfile:
            cds_with_absent_proteins_id_list.write(cds_with_no_pep_idfile)

    def get_transcript_to_pep_accordance_from_files(self, transcript_file, pep_file, output_file, transcript_with_no_pep_idfile=None, format="fasta",
                                                    verbose=True, parsing_mode="parse", index_file_suffix="tmp.idx",
                                                    transcript_with_several_proteins_idfile=None,
                                                    genetic_code_table=1, include_id_check=False):

        transcript_dict = self.parse_seq_file(transcript_file, mode=parsing_mode, format=format,
                                              index_file="cds_%s" % index_file_suffix)
        pep_dict = self.parse_seq_file(pep_file, mode=parsing_mode, format=format,
                                       index_file="pep_%s" % index_file_suffix)

        transcript_pep_accordance_dict, transcript_with_absent_proteins_id_list, transcript_with_several_proteins = \
            self.get_transcript_to_pep_accordance(transcript_dict, pep_dict, verbose=verbose,
                                                  parsing_mode=parsing_mode,
                                                  genetic_code_table=genetic_code_table,
                                                  include_id_check=include_id_check)
        transcript_pep_accordance_dict.write(output_file, splited_values=True)
        if transcript_with_no_pep_idfile:
            transcript_with_absent_proteins_id_list.write(transcript_with_no_pep_idfile)

        if transcript_with_several_proteins_idfile:
            transcript_with_several_proteins.write(transcript_with_several_proteins_idfile)

    @staticmethod
    def extract_exon_lengths(record_dict):
        output_list = []
        for record_id in record_dict:
            for feature in record_dict[record_id].features:
                taxonomy = ";".join(record_dict[record_id].annotations["taxonomy"])
                species = record_dict[record_id].annotations["organism"]
                if feature.type == "mRNA" or feature.type == "transcript":
                    product = ";".join(feature.qualifiers["product"]) if "product" in feature.qualifiers else "."
                    transcript_id = ",".join(feature.qualifiers["transcript_id"]) if "transcript_id" in feature.qualifiers else "."
                    strand = feature.location.strand
                    exon_lengths = []
                    #print feature.sub_features
                    #print feature.location
                    #print feature.location.start
                    for location in feature.location.parts:
                        #print location
                        exon_len = location.end - location.start
                        exon_lengths.append(exon_len)

                    output_list.append([species, taxonomy, record_id, transcript_id, product, strand, exon_lengths])
        return output_list

    def extract_exon_lengths_from_genbank_file(self, input_file, output_file):
        record_dict = SeqIO.index_db("tmp.idx", input_file, format="genbank")
        data_list = self.extract_exon_lengths(record_dict)

        with open(output_file, "w") as out_fd:
            out_fd.write("#species\ttaxonomy\trecord_id\ttranscript_id\tproduct\tstrand\texon_length\n")
            for entry in data_list:
                out_fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entry[0], entry[1], entry[2], entry[3], entry[4],
                                                               str(entry[5]), ",".join(map(str, entry[6]))))

        os.remove("tmp.idx")

    @staticmethod
    def get_degenerate_codon_set(genetic_code_table, type="dna"):
        if type == "dna":
            nucleotides = ["A", "T", "G", "C"]
        elif type == "rna":
            nucleotides = ["A", "U", "G", "C"]
        degenerate_codon_set = set()
        for pos_1 in nucleotides:
            for pos_2 in nucleotides:
                codon = pos_1 + pos_2 + "N"
                if Seq(codon).translate(table=genetic_code_table) != "X":
                    degenerate_codon_set.add(codon)
        return degenerate_codon_set

    @staticmethod
    def extract_introns_from_transcripts(record_dict, transcript_id_white_list=None):
        intron_dict = OrderedDict()
        unknown_transcript_index = 1
        for record_id in record_dict:
            for feature in record_dict[record_id].features:
                used_transcript_id = None

                if feature.type == "mRNA" or feature.type == "transcript":
                    if transcript_id_white_list:
                        for transcript_id in feature.qualifiers["transcript_id"]:
                            if transcript_id in transcript_id_white_list:
                                used_transcript_id = transcript_id
                                break

                        if used_transcript_id is None:
                            continue
                    #product = ";".join(feature.qualifiers["product"]) if "product" in feature.qualifiers else "."
                    transcript_id = used_transcript_id if used_transcript_id else feature.qualifiers["transcript_id"][0] if "transcript_id" in feature.qualifiers else None
                    if transcript_id is None:
                        transcript_id = "unkown_transcript_%i" % unknown_transcript_index
                        unknown_transcript_index += 1

                    strand = feature.location.strand

                    number_of_exons = len(feature.location.parts)
                    for i in range(0, number_of_exons - 1):
                        intron_location = FeatureLocation(feature.location.parts[i].end,
                                                          feature.location.parts[i+1].start,
                                                          strand=strand)
                        if strand >= 0:
                            previous_exon_number = i + 1
                            following_exon_number = i + 2
                            previous_exon_len = len(feature.location.parts[i])
                            following_exon_len = len(feature.location.parts[i+1])
                        else:
                            previous_exon_number = number_of_exons - i - 1
                            following_exon_number = number_of_exons - i
                            previous_exon_len = len(feature.location.parts[i+1])
                            following_exon_len = len(feature.location.parts[i])

                        intron_id = "%s_intron_%i-%i_between_exons_%i-%i" % (transcript_id, previous_exon_number,
                                                                             following_exon_number, previous_exon_len,
                                                                             following_exon_len)
                        intron_record = SeqRecord(intron_location.extract(record_dict[record_id].seq), id=intron_id,
                                                  description=record_dict[record_id].annotations["organism"])

                        intron_dict[intron_id] = intron_record

        return intron_dict

    def extract_introns_from_transcripts_from_genbank_files(self, list_of_genbank_files, output_file,
                                                            transcript_id_white_list=None):
        record_dict = SeqIO.index_db("tmp.idx", list_of_genbank_files, format="genbank")
        intron_dict = self.extract_introns_from_transcripts(record_dict,
                                                            transcript_id_white_list=transcript_id_white_list)
        SeqIO.write(self.record_from_dict_generator(intron_dict), output_file, format="fasta")

        os.remove("tmp.idx")

    @staticmethod
    def get_general_statistics(input_file, file_format, output_file="statistics.t", write_to_file=False):
        record_dict = SeqIO.to_dict(SeqIO.parse(input_file, file_format))
        statistics_dict = {}
        for record_id in record_dict:
            statistics_dict[record_id] = [len(record_dict[record_id].seq), len(record_dict[record_id].features)]
        if write_to_file:
            fd = open(output_file, "w")
            metadata = "#Totaly records\t%i\n" % len(record_dict)
            fd.write(metadata)
            for record_id in sorted(list(record_dict.keys())):
                fd.write("%s\t%i\t%i\n" % (record_id, statistics_dict[record_id][0], statistics_dict[record_id][1]))
            fd.close()
        return statistics_dict

    def get_multifile_general_statistics(self, input_file_list, file_format, output_file="statistics.t", 
                                         output_file2="file_statistics.t", names_list=None, write_to_file=False):
        statistics_dict = {}
        names_dict = {}
        if names_list:
            for i in range(0, len(input_file_list)):
                names_dict[input_file_list[i]] = names_list[i]
        else:
            for i in range(0, len(input_file_list)):
                names_dict[input_file_list[i]] = input_file_list[i]
    
        for file_entry in input_file_list:
            print("Gathering statistics from %s" % file_entry)
            statistics_dict[names_dict[file_entry]] = self.get_general_statistics(file_entry, file_format)
    
        if write_to_file:
            fd_file = open(output_file2, "w")
            fd = open(output_file, "w")
            metadata = "#Totaly files\t%i\n" % len(input_file_list)
            fd.write(metadata)
            fd_file.write(metadata)
            for entry in sorted(list(statistics_dict.keys())):
                fd.write("%s\t%i\n" % (entry, len(statistics_dict[entry])))
                fd_file.write("%s\t%i\n" % (entry, len(statistics_dict[entry])))
                for record_id in sorted(list(statistics_dict[entry].keys())):
                    fd.write("\t%s\t%i\t%i\n" % (record_id, statistics_dict[entry][record_id][0], statistics_dict[entry][record_id][1]))
            fd.close()
            fd_file.close()
        return statistics_dict
    
    def get_statistics(self, filelist, index_filename, taxonomy_file_prefix="taxonomy", filetype="genbank"):
        record_dict = SeqIO.index_db(index_filename, filelist, filetype)
        number_of_records = len(record_dict)
        print("Gathering statistics...")
        print("Total records:\t %i" % number_of_records)
        taxonomy_dict = {}
        # constructing nested taxonomy summary as nested dictionary
        print("Constructing taxonomy distribution...")
    
        for record_id in record_dict:
            temp_dict = taxonomy_dict

            for taxon in record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]:
                if taxon not in temp_dict:
                    temp_dict[taxon] = {"Total": 1}
                else:
                    temp_dict[taxon]["Total"] += 1
    
                temp_dict = temp_dict[taxon]

        with open(taxonomy_file_prefix + ".pickle", 'wb') as fd:
            pickle.dump(taxonomy_dict, fd)
        self.print_formatted_nested_dict(taxonomy_dict, taxonomy_file_prefix + ".tax", indent_type="\t")
    
        return record_dict, taxonomy_dict
    
    @staticmethod
    def get_taxonomy(record_dict, output_prefix="taxonomy"):
        number_of_records = len(record_dict)
        print("Gathering taxonomy statistics...")
        print("Total records:\t %i" % number_of_records)
        genus_dict = {}
        family_dict = {}
        order_dict = {}
        # constructing nested taxonomy summary as nested dictionary
        print("Constructing taxonomy distribution...")
    
        for record_id in record_dict:
            if record_dict[record_id].annotations["taxonomy"][-1] not in genus_dict:
                genus_dict[record_dict[record_id].annotations["taxonomy"][-1]] = [record_id]
            else:
                genus_dict[record_dict[record_id].annotations["taxonomy"][-1]].append(record_id)
    
            for taxon in record_dict[record_id].annotations["taxonomy"]:
                if taxon[-6:] == "formes":
                    if taxon not in order_dict:
                        order_dict[taxon] = [record_id]
                    else:
                        order_dict[taxon].append(record_id)
                if taxon[-4:] == "idae":
                    if taxon not in family_dict:
                        family_dict[taxon] = [record_id]
                    else:
                        family_dict[taxon].append(record_id)
    
        def write_taxa_dict(taxa_dict, taxa_file):
            with open(taxa_file, "w") as fd:
                for taxon in taxa_dict:
                    fd.write("%s\t%i\t" % (taxon, len(taxa_dict[taxon])) + ",".join(taxa_dict[taxon]) + "\n")
    
        write_taxa_dict(genus_dict, output_prefix + "_genus.tax")
        write_taxa_dict(family_dict, output_prefix + "_family.tax")
        write_taxa_dict(order_dict, output_prefix + "_order.tax")
    
        return genus_dict, family_dict, order_dict
    
    def get_taxonomy_from_genbank_files(self, filelist, index_filename):
        record_dict = SeqIO.index_db(index_filename, filelist, "genbank")
        return self.get_taxonomy(record_dict)
    
    @staticmethod
    def count_species(record_dict, output_filename="count_species.count"):
        species_count_dict = {}
        i = 1
        print("Counting species...")
        for record_id in record_dict:
            organism = record_dict[record_id].annotations['organism']
            if organism not in species_count_dict:
                species_count_dict[organism] = [1, [record_id]]
            else:
                species_count_dict[organism][0] += 1
                species_count_dict[organism][1].append(record_id)
            i += 1
        number_of_species = len(species_count_dict)
        print ("Total %i species were found" % number_of_species)
        fd = open(output_filename, "w")
        fd.write("#Number of species\t%i\n" % number_of_species)
        for organism in species_count_dict:
            fd.write(organism + "\t%i\t%s\n" % (species_count_dict[organism][0], 
                                                ",".join(species_count_dict[organism][1])))
        return number_of_species, species_count_dict

    def count_species_from_file(self, sequence_file, format="genbank", output_filename="count_species.count"):
        record_dict = SeqIO.index_db("tmp.idx", sequence_file, format=format)
        self.count_species(record_dict, output_filename=output_filename)
        os.remove("tmp.idx")

    @staticmethod
    def get_id_to_species_accordance(record_dict, output="id_to_species.accordance"):
        accordance_dict = {}
        for record_id in record_dict:
            organism = record_dict[record_id].annotations['organism'] if "organism" in record_dict[record_id].annotations else "."
            accordance_dict[record_id] = organism

        with open(output, "w") as out_fd:
            for record_id in accordance_dict:
                print(record_id, accordance_dict[record_id])
                out_fd.write("%s\t%s\n" % (record_id, accordance_dict[record_id]))

        return accordance_dict

    def get_id_to_species_accordance_from_file(self, sequence_file, format="genbank",
                                               output="id_to_species.accordance"):
        record_dict = SeqIO.index_db("tmp.idx", sequence_file, format=format)
        self.get_id_to_species_accordance(record_dict, output=output)
        os.remove("tmp.idx")

    @staticmethod
    def split_records_by_taxa_level(record_dict, prefix, taxa_level=1, filetype="genbank"):
        """taxa levels starts from 0"""
        taxa_dict = {}
        print("Splitting records by %i taxa level" % taxa_level)
        for record_id in record_dict:
            if len(record_dict[record_id].annotations["taxonomy"]) > taxa_level:
                taxon = record_dict[record_id].annotations["taxonomy"][taxa_level]
                if taxon not in taxa_dict:
                    taxa_dict[taxon] = [record_id]
                else:
                    taxa_dict[taxon].append(record_id)
            else:
                print(record_id,
                      "Not classified for selected taxa level",
                      "Taxonomy:",
                      record_dict[record_id].annotations["taxonomy"],
                      record_dict[record_id].annotations["organism"])
    
        def taxon_generator():
            for record_id in taxa_dict[taxon]:
                yield record_dict[record_id]
    
        for taxon in taxa_dict:
            taxon_name = taxon.replace(" ", "_")
            filename = str(prefix) + "_" + str(taxa_level) + "_" + str(taxon_name) + ".gb"
            print(filename, len(taxa_dict[taxon]))
            SeqIO.write(taxon_generator(), filename, filetype)
    
        number_of_taxa = len(taxa_dict)
        print("Totaly %i taxa of %i level were found" % (number_of_taxa, taxa_level))
        return number_of_taxa
    
    @staticmethod
    def print_formatted_nested_dict(nested_dict, output_file, indent_type="\t"):
        #at moment only \t indent is supported
    
        def fwrite_nested_dict(nested_dict, fd, indent_type, indent_counter):
            for taxon in nested_dict:
                if taxon != "Total":
                    #print(nested_dict[taxon])
                    fd.write(indent_type*indent_counter + taxon + "\t%i\n" % nested_dict[taxon]["Total"])
                    fwrite_nested_dict(nested_dict[taxon], fd, indent_type, indent_counter + 1)
    
        fd = open(output_file, "w")
        fwrite_nested_dict(nested_dict, fd, indent_type, 0)
        fd.close()
        
    def filter_by_taxa(self, record_dict,
                       taxa_list,
                       filter_type="white_list",
                       output_filename="filtered_by_taxa.gb",
                       store_filtered_out=True,
                       filtered_out_filename="filtered_out_by_taxa.gb",
                       output_type="genbank",
                       return_record_generator=False):
        print("Filtering by taxa...")
        print("Totaly %i records" % len(record_dict))
    
        filtered_id_list = []
        filtered_out_id_list = []
        taxa_set = set(taxa_list)
        if filter_type == "white_list":
            print("Taxa to be retained: %s" % ", ".join(taxa_list))
            for record_id in record_dict:
                #print(record_id)
                if set(record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]) & taxa_set:
                    filtered_id_list.append(record_id)
                else:
                    filtered_out_id_list.append(record_id)
        elif filter_type == "black_list":
            print("Taxa to be filtered out: %s" % ", ".join(taxa_list))
            for record_id in record_dict:
                for filter_entry in taxa_set:
                    for taxa_entry in record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]:
                        if filter_entry in taxa_entry:
                            filtered_out_id_list.append(record_id)
                            break
                    else:
                        continue
                    break
                else:
                    filtered_id_list.append(record_id)
        SeqIO.write(self.record_by_id_generator(record_dict, filtered_id_list), output_filename, output_type)
        if store_filtered_out:
            SeqIO.write(self.record_by_id_generator(record_dict, filtered_out_id_list), filtered_out_filename, output_type)
        number_of_retained = len(filtered_id_list)
        print("Retained %i records" % number_of_retained)
    
        if return_record_generator:
            return self.record_by_id_generator(record_dict, filtered_id_list)
    
    def split_by_taxa(self, record_dict, taxa_list, output_suffix, output_type="genbank", outfiltered_name="outfiltered"):
        print("Spliting taxa...")
        taxa_dict = {}
        for taxa in taxa_list:
            taxa_dict[taxa] = []
        taxa_dict[outfiltered_name] = []
    
        for record_id in record_dict:
            record_taxa = record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]
            for taxa in taxa_list:
                if taxa in record_taxa:
                    taxa_dict[taxa].append(record_id)
                    break
            else:
                taxa_dict[outfiltered_name].append(record_id)
    
        if output_type == "genbank":
            out_extension = ".gb"
    
        for taxa in taxa_dict:
            print("Writing %s taxa" % taxa)
            SeqIO.write(self.record_by_id_generator(record_dict, taxa_dict[taxa]),
                        taxa + output_suffix + out_extension, output_type)
        return taxa_dict

    def filter_by_source(self, record_dict,
                         source_list,
                         filter_type="white_list",
                         output_filename="filtered_by_source.gb",
                         output_type="genbank",
                         return_record_generator=False):
        print("Filtering by source")
        print("Totaly %i records" % len(record_dict))
        filtered_id_list = []
        filtered_out_id_list = []
        source_set = set(source_list)
        if filter_type == "white_list":
            print("Sources to be retained: %s" % ", ".join(source_list))
            for record_id in record_dict:
                if ("source" not in record_dict[record_id].annotations) or (not (record_dict[record_id].annotations["source"])):
                    continue
                """
                products = ""
                for feature in record_dict[record_id].features:
                    if "product" in feature.qualifiers:
                        products += " " + feature.qualifiers["product"]
                """
                if set(record_dict[record_id].annotations["source"].split()) & source_set: # or (product_set & source_set):
                    filtered_id_list.append(record_id)
                else:
                    filtered_out_id_list.append(record_id)
        elif filter_type == "black_list":
            print("Sources to be filtered out: %s" % ", ".join(source_list))
            for record_id in record_dict:
                if "source" not in record_dict[record_id].annotations:
                    continue
                """
                product_set = set([])
                for feature in record_dict[record_id].features:
                    if "product" in feature.qualifiers:
                        product_set = set(feature.qualifiers["product"])
                """
                if not (set(record_dict[record_id].annotations["source"].split()) & source_set):# or (product_set & source_set)):
                    filtered_id_list.append(record_id)
                else:
                    filtered_out_id_list.append(record_id)
        SeqIO.write(self.record_by_id_generator(record_dict, filtered_id_list), output_filename, output_type)
        SeqIO.write(self.record_by_id_generator(record_dict, filtered_out_id_list), "filtered_out_by_source.gb", output_type)
        number_of_retained = len(filtered_id_list)
        print("Retained %i records" % number_of_retained)
    
        if return_record_generator:
            return self.record_by_id_generator(record_dict, filtered_id_list)
    
    @staticmethod
    def parse_counts_file(species_count_file):
        fd = open(species_count_file, "r")
        number_of_species = fd.readline().strip().split("\t")[-1]
        species_count_dict = {}
        for line in fd:
            line_list = line.strip().split("\t")
            species_count_dict[line_list[0]] = [int(line_list[1]), line_list[2].split(",")]
        fd.close()
        return number_of_species, species_count_dict

    def sort_by_number(self, species_count_file, output_filename_prefix, counts_list=[1, 5]):
        number_of_species, species_count_dict = self.parse_counts_file(species_count_file)
        length_of_counts = len(counts_list)
        fd_list = [open(output_filename_prefix + "_%i_%i_species.count" % (counts_list[i], counts_list[i+1] - 1), "w") for i in range(0, length_of_counts-1)]
        fd_list.append(open(output_filename_prefix + "_%i+_species.count" % (counts_list[-1]), "w"))
        number_of_species = [0 for i in range(0, length_of_counts)]
        species_record_list = [[] for i in range(0, length_of_counts)]
    
        for species in species_count_dict:
            for i in range(0, length_of_counts-1):
                if counts_list[i] <= species_count_dict[species][0] < counts_list[i+1]:
                    number_of_species[i] += 1
                    species_record_list[i].append(species)
                    break
            else:
                number_of_species[-1] += 1
                species_record_list[-1].append(species)

        for i in range(0, length_of_counts):
            fd_list[i].write("#Number of species\t%i\n" % number_of_species[i])
            for species in species_record_list[i]:
                fd_list[i].write(species + "\t%i\t%s\n" % (species_count_dict[species][0], ",".join(species_count_dict[species][1])))
            fd_list[i].close()

    @staticmethod
    def renamed_records_generator(record_dict, syn_dict=None, expression=None, clear_description=False,
                                  store_old_name_in_description=False):
        for record_id in record_dict:
            if syn_dict:
                if record_id not in syn_dict:
                    print("%s was not renamed" % record_id)
                    yield record_dict[record_id]
                    continue
                else:
                    record = deepcopy(record_dict[record_id])
                    record.id = syn_dict[record_id]

            elif expression:
                record = deepcopy(record_dict[record_id])
                record.id = expression(record_id)

            if clear_description:
                record.description = ""
            if store_old_name_in_description:
                record.description += " old_name=%s" % record_id

            yield record

    def rename_records_from_files(self, input_file, output_file, synonyms_file=None, format="fasta", header=False,
                                  separator="\t", key_index=0, value_index=1, syn_expression=None, comments_prefix=None,
                                  clear_description=False, record_id_expression=None, store_old_name_in_description=False,
                                  parsing_mode="parse"):
        if synonyms_file and record_id_expression:
            raise ValueError("Both synonyms file and record id expression were set")
        elif (not synonyms_file) and (not record_id_expression):
            raise ValueError("Neither synonyms file nor record id expression were set")
        elif synonyms_file:
            syn_dict = SynDict()
            syn_dict.read(synonyms_file, header=header, separator=separator, key_index=key_index,
                          value_index=value_index, expression=syn_expression, comments_prefix=comments_prefix)
        else:
            syn_dict = None

        record_dict = self.parse_seq_file(input_file, parsing_mode, format=format, index_file="temp.idx")
        #SeqIO.index_db("temp.idx", input_file, format=format)

        SeqIO.write(self.renamed_records_generator(record_dict, syn_dict=syn_dict, expression=record_id_expression,
                                                   store_old_name_in_description=store_old_name_in_description,
                                                   clear_description=clear_description),
                    output_file, format=format)
        if parsing_mode == "index_db":
            os.remove("temp.idx")

    def rename_records_by_sequential_ids_from_files(self, input_file, output_file, output_syn_file, format="fasta",
                                                    clear_description=False, record_id_prefix="SEQ",
                                                    length_of_numerical_part=8, parse_mode="parse",
                                                    index_file="temp.idx"):

        record_dict = self.parse_seq_file(input_file, mode=parse_mode, format=format, index_file=index_file)

        syn_dict = SynDict()
        id_template = "%s%%0%ii" % (record_id_prefix, length_of_numerical_part)
        i = 1
        for record_id in record_dict:
            syn_dict[record_id] = id_template % i
            i += 1
        syn_dict.write(output_syn_file)

        SeqIO.write(self.renamed_records_generator(record_dict, syn_dict=syn_dict, expression=None,
                                                   clear_description=clear_description),
                    output_file, format=format)

    @staticmethod
    def trim_cds_and_remove_terminal_stop_codons_generator(record_dict, stop_codons_list=("TGA", "TAA", "TAG")):

        stop_codons = []
        for stop_codon in stop_codons_list:
            stop_codons.append(stop_codon.upper().replace("U", "T"))

        for record_id in record_dict:
            new_record = deepcopy(record_dict[record_id])
            seq_length = len(new_record.seq)
            remainder = seq_length % 3
            if remainder == 0:
                if str(new_record.seq[-3:]).upper() in stop_codons_list:
                    trim = 3
                else:
                    trim = 0
            else:
                if str(new_record.seq[-3+remainder:remainder]).upper() in stop_codons:
                    trim = 3 + remainder
                else:
                    trim = remainder
            if trim > 0:
                new_record.seq = new_record.seq[:-trim]

            yield new_record

    def trim_cds_and_remove_terminal_stop_codons(self, input_file, output_file, stop_codons_list=("TGA", "TAA", "TAG")):

        record_dict = SeqIO.index_db("tmp.idx", input_file, format="fasta")

        SeqIO.write(self.trim_cds_and_remove_terminal_stop_codons_generator(record_dict,
                                                                            stop_codons_list=stop_codons_list),
                    output_file, format="fasta")
        os.remove("tmp.idx")

    # --------------------Search--------------------------
    #              !!!!IMPORTANT!!!!!!
    # in this section python notation for coordinates inside sequence is used
    @staticmethod
    def find_homopolymer_end(seq, nucleotide, seq_length, start, search_type="perfect",
                           max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2):
        shift_list = []
        number_of_insertions = 0
        total_insertion_length = 0
        insertion_length = 0
        i = start + 1
        if search_type == 'perfect':
            while i < seq_length:
                if seq[i] != nucleotide:
                    return i, None
                i += 1
            return seq_length, None
        else:
            while i < seq_length:
                if seq[i] != nucleotide:
                    if seq[i-1] == nucleotide:
                        shift_list.append(i)
                        insertion_length = 1
                    else:
                        insertion_length += 1
                    number_of_insertions += 1
                    total_insertion_length += 1
                    if number_of_insertions > max_number_of_insertions or insertion_length > max_single_insert_size:
                        end = shift_list[-1]
                        break
                    if max_total_insert_length:
                        if total_insertion_length > max_total_insert_length:
                            end = shift_list[-1]
                            break
                i += 1
            else:
                if seq[-1] == nucleotide:
                    end = seq_length
                else:
                    end = shift_list[0]

            new_start = shift_list[0] if shift_list else end
            return end, new_start

    @staticmethod
    def make_region_bed_file(record_dict, output_file, white_list=None, black_list=None,
                             output_format="0-based", min_len=None, max_len=None):
        with open(output_file, "w") as out_fd:
            for record_id in record_dict:
                if white_list:
                    if record_id not in white_list:
                        continue
                if black_list:
                    if record_id in black_list:
                        continue
                if min_len:
                    if len(record_dict[record_id]) < min_len:
                        continue
                if max_len:
                    if len(record_dict[record_id]) > max_len:
                        continue
                bed_string = "%s" % record_id
                if output_format == "0-based":
                    bed_string += "\t0\t%i\n" % (len(record_dict[record_id]) - 1)
                elif output_format == "1-based":
                    bed_string += "\t1\t%i\n" % len(record_dict[record_id])
                else:
                    raise ValueError("Unrecognized bed format %s" % output_format)
                out_fd.write(bed_string)

    def make_region_bed_file_from_file(self, seq_file, output_file, white_id_file=None, black_id_file=None,
                                       output_format="0-based", input_format="fasta", min_len=None, max_len=None,
                                       parsing_mode="index_db", index_file="tmp.idx", retain_index=False):
        record_dict = self.parse_seq_file(seq_file, parsing_mode, format=input_format, index_file=index_file)#SeqIO.index_db("tmp.idx", seq_file, format=input_format)
        black_id_list = IdList()
        white_id_list = IdList()
        if black_id_file:
            black_id_list.read(black_id_file)
        if white_id_file:
            white_id_list.read(white_id_file)

        self.make_region_bed_file(record_dict, output_file, white_list=white_id_list, black_list=black_id_list,
                                  output_format=output_format, min_len=min_len, max_len=max_len)
        if (parsing_mode == "index_db") and (not retain_index):
            os.remove(index_file)

    def find_homopolymers(self, seq, nucleotide, min_size=5, search_type="perfect",
                          max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2):
        # search types:
        #   perfect - search only for perfect homopolymers, all options other than min_size are ignored
        #   non_perfect - search for non_perfect homopolymers with max_single_insert_size, max_total_insert_length
        #                 and max_number_of_insertions
        seq_length = len(seq)
        i = 0
        homopolymers_coords = []
        homopolymers_lengthes = []
        #print(seq_length)
        prev_end = 0
        while i < seq_length:
            if seq[i] == nucleotide:
                end, new_start = self.find_homopolymer_end(seq, nucleotide, seq_length, i, search_type=search_type,
                                                           max_single_insert_size=max_single_insert_size,
                                                           max_total_insert_length=max_total_insert_length,
                                                           max_number_of_insertions=max_number_of_insertions)
                # print(end, i)
                length = end - i
                if homopolymers_coords:
                    prev_end = homopolymers_coords[-1][1]

                if length >= min_size and end > prev_end:
                    homopolymers_coords.append((i, end))
                    homopolymers_lengthes.append(length)

                i = end
                if new_start is not None:
                    i = new_start
                continue
            i += 1

        homopolymers_lengthes = np.array(homopolymers_lengthes)
        return homopolymers_coords, homopolymers_lengthes

    @staticmethod
    def calculate_assembly_stats(record_dict, thresholds_list=(0, 500, 1000), seq_len_file=None):
        Ns_number = 0
        for contig_id in record_dict:
            Ns_number += record_dict[contig_id].seq.count("N") + record_dict[contig_id].seq.count("n")

        # for record in record_dict:
        #    print "%s\t%i" % (record, len(record_dict[record].seq))

        length_array = np.array(sorted([len(record_dict[record].seq) for record in record_dict], reverse=True))
        if seq_len_file:
            np.savetxt(seq_len_file, length_array, fmt="%i")
        total_length = sum(length_array)
        longest_contig = length_array[0]

        right_bin = int(math.log10(longest_contig)) + 1
        bins = [10**i for i in range(0, right_bin + 1)]
        bins = bins
        contig_cumulative_length_values = [0 for i in range(0, right_bin)]
        contig_number_values = [0 for i in range(0, right_bin)]

        #print length_array
        for contig_len in length_array:
            #print contig_len
            len_power = int(math.log10(contig_len))
            contig_cumulative_length_values[len_power] += contig_len
            contig_number_values[len_power] += 1

        #print(bins)
        #print(contig_cumulative_length_values)
        #print(contig_number_values)
        L50_dict = OrderedDict()
        N50_dict = OrderedDict()
        length_dict = OrderedDict()
        # TODO: make calculations of all N50s in one run
        for threshold in thresholds_list:
            length_above_threshold = 0
            N50 = 0
            for length in length_array:
                if length >= threshold:
                    length_above_threshold += length
                else:
                    break
            half_length = int(length_above_threshold/2)
            half_length = half_length if half_length*2 == length_above_threshold else half_length + 1

            tmp = 0
            number_of_contigs = 0
            tmp_length = 0
            for length in length_array:
                """
                if tmp < half_length:
                    number_of_contigs += 1
                    tmp += length
                else:
                    N50 = length
                    break
                """
                number_of_contigs += 1
                tmp += length
                if tmp >= half_length:
                    N50 = length
                    break
            N50_dict[threshold] = N50
            L50_dict[threshold] = number_of_contigs
            for length in length_array:
                if length >= threshold:
                    tmp_length += length
                else:
                    break
            length_dict[threshold] = tmp_length

        return length_array, N50_dict, L50_dict, length_dict, total_length, longest_contig, Ns_number, bins, contig_cumulative_length_values, contig_number_values

    def get_random_species_genomes(self,
                                   record_dict,
                                   output_file,
                                   selected_species_file="selected_species.t",
                                   count_species_file=None,
                                   output_type="fasta",
                                   prev_id_dict={}):
        print("Extracting random genomes(one per species)...")
        #return ids of random genomes, one per species

        if not count_species_file:
            self.count_species(record_dict, output_filename="count_species.count")
            fd = open("count_species.count", "r")
        else:
            fd = open(count_species_file, "r")
        fd_species = open(selected_species_file, "w")
        fd_species.write("#species\tid\n")
        number_of_species = int(fd.readline().strip().split("\t")[1])
        print("Totaly %s species " % number_of_species)
        random_ids = []
        retained_ids = []
        for line in fd:
            line_list = line.strip().split("\t")
            species = line_list[0]
            if species in prev_id_dict:
                if prev_id_dict[species] in record_dict:
                    new_id = prev_id_dict[species]
                    retained_ids.append(prev_id_dict[species])
            else:
                num_of_genomes = int(line_list[1])
                id_list = line_list[2].split(",")
                random_number = randint(1, num_of_genomes)
                new_id = id_list[random_number-1]
            random_ids.append(new_id)
            fd_species.write("%s\t%s\n" % (species, new_id))
        fd.close()
        fd_species.close()
        random_ids = set(random_ids)
        SeqIO.write([record_dict[id_entry] for id_entry in random_ids], open(output_file, "w"), output_type)
        #print(len(retained_ids), retained_ids)
        return random_ids

    def get_random_species_genomes_from_genbank_file(self, input_gb_files, output_prefix, output_type="genbank"):

        record_dict = SeqIO.index_db("tmp.idx", self.make_list_of_path_to_files(input_gb_files),
                                     format="genbank")
        #count_species_file = "%s.species_counts" % output_prefix
        output_file = "%s.extracted.%s" % (output_prefix, output_type)
        self.get_random_species_genomes(record_dict, output_file,
                                        output_type=output_type, prev_id_dict={})

        os.remove("tmp.idx")

    @staticmethod
    def get_protein_segments_from_CDS_location(CDS_location, number_of_end_nucleotides_to_ignore=0):
        cds_segments_len_list = [len(segment_location) for segment_location in CDS_location.parts]
        #print("CDS segment list")
        #print cds_segments_len_list
        cds_segments_len_list[-1] -= number_of_end_nucleotides_to_ignore
        #print cds_segments_len_list
        tmp_len = 0
        next_segment_start = 0
        left_phase = 0
        right_phase = 0

        protein_segments_list = []
        protein_segment_phases_list = []
        for i in range(0, len(cds_segments_len_list)):
            tmp_len += cds_segments_len_list[i]
            full_triplet_cummulative_number = int(tmp_len / 3)
            remainder = tmp_len % 3

            right_phase = remainder
            segment_location = FeatureLocation(start=next_segment_start, end=full_triplet_cummulative_number)
            segment_phase = (left_phase, right_phase)
            if remainder == 0:
                next_segment_start = full_triplet_cummulative_number
                left_phase = 0
            else:
                next_segment_start = full_triplet_cummulative_number + 1
                left_phase = 3 - remainder

            protein_segments_list.append(segment_location)
            protein_segment_phases_list.append(segment_phase)
        #python-based coordinates
        return protein_segments_list, protein_segment_phases_list

    @staticmethod
    def make_string_from_nested_list_of_ints(input_list, first_lvl_separator=",", second_lvl_separator=";"):
        return second_lvl_separator.join([first_lvl_separator.join(map(str, coordinates)) for coordinates in input_list])

    def get_protein_marking_by_exons_from_genbank(self, record_dict,
                                                  output_prefix,
                                                  protein_id_field_in_cds_feature="protein_id",
                                                  gene_id_field_in_cds_feature="gene",
                                                  translation_field_in_cds_feature="translation",
                                                  genetic_code_table=1, stop_codon_symbol="*", verbose=True):
        """
        def extract_feature_seq(feature, record):
            seq = ""
            for part in feature.location.parts:
                tmp = record.seq[part.start:part.end]
                tmp = tmp.reverse_complement() if feature.location.strand == -1 else tmp
                seq += str(tmp)
            return seq
        """

        output_full_header_list = ["protein_id",
                                   "gene_id",
                                   "protein_length",
                                   "CDS_exon_number",
                                   "CDS_exon_coordinates",
                                   "CDS_strand",
                                   "protein_segments",
                                   "protein_segment_phases",
                                   "CDS_seq",
                                   "protein_seq"]

        output_simple_header_list = ["protein_id",
                                     "protein_length",
                                     "segment_number",
                                     "protein_segments",
                                     "protein_segment_phases"
                                     ]

        output_full_header = "\t".join(output_full_header_list) + "\n"
        output_simple_header = "\t".join(output_simple_header_list) + "\n"

        output_full_file = "%s.full" % output_prefix
        output_simple_file = "%s.simple" % output_prefix

        out_full_fd = open(output_full_file, "w")
        out_full_fd.write(output_full_header)

        out_simple_fd = open(output_simple_file, "w")
        out_simple_fd.write(output_simple_header)

        ignored_protein_counter = 0
        for record_id in record_dict:
            if verbose:
                print("Handling %s" % record_id)
            for feature in record_dict[record_id].features:

                if feature.type == "CDS":
                    protein_id = feature.qualifiers[protein_id_field_in_cds_feature][0] if protein_id_field_in_cds_feature in feature.qualifiers else None
                    gene_id = feature.qualifiers[gene_id_field_in_cds_feature][0] if gene_id_field_in_cds_feature in feature.qualifiers else None
                    translation = feature.qualifiers[translation_field_in_cds_feature][0] if translation_field_in_cds_feature in feature.qualifiers else None

                    translation = translation[:-1] if translation[-1] == stop_codon_symbol else translation

                    if verbose:
                        print("\t" + protein_id)
                    #print gene_id
                    #print translation
                    CDS_len = len(feature.location)
                    CDS = feature.extract(record_dict[record_id]).seq
                    CDS_translation = CDS.translate(to_stop=True, table=genetic_code_table,
                                                    stop_symbol=stop_codon_symbol)
                    #print CDS_translation

                    protein_length = len(translation)

                    #print feature.location
                    #print feature.location.strand
                    #print feature.location.parts
                    segments_lengths = [len(segment_location) for segment_location in feature.location.parts]
                    #print segments_lengths
                    CDS_exon_number = len(feature.location.parts)
                    number_of_end_nucleotides_to_ignore = CDS_len - protein_length*3
                    protein_segments, protein_segment_phases_list = self.get_protein_segments_from_CDS_location(feature.location,
                                                                                                                number_of_end_nucleotides_to_ignore=number_of_end_nucleotides_to_ignore)
                    #print protein_segments
                    #print CDS_len
                    #print float(len(feature.location))/3.0
                    #print len(translation) if translation else None
                    #print "\n\n"

                    CDS_exon_coordinates = [[location.start, location.end] for location in feature.location.parts]
                    CDS_exon_coordinates[-1][-1] -= number_of_end_nucleotides_to_ignore
                    CDS_exon_coordinates_for_string = deepcopy(CDS_exon_coordinates)
                    for i in range(0, len(CDS_exon_coordinates)):
                        CDS_exon_coordinates_for_string[i][0] += 1
                    CDS_exon_coordinates_str = self.make_string_from_nested_list_of_ints(CDS_exon_coordinates_for_string)

                    protein_segments_for_string = deepcopy(protein_segments)
                    for i in range(0, len(protein_segments_for_string)):
                        protein_segments_for_string[i] = [protein_segments_for_string[i].start + 1, protein_segments_for_string[i].end]
                    output_full_string = "%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%s\t%s\n" % (protein_id,
                                                                                       gene_id,
                                                                                       protein_length,
                                                                                       CDS_exon_number,
                                                                                       CDS_exon_coordinates_str,
                                                                                       "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else ".",
                                                                                       self.make_string_from_nested_list_of_ints(protein_segments_for_string),
                                                                                       self.make_string_from_nested_list_of_ints(protein_segment_phases_list),
                                                                                       str(CDS)[:-number_of_end_nucleotides_to_ignore],
                                                                                       str(CDS_translation))

                    output_simple_string = "%s\t%i\t%i\t%s\t%s\n" % (protein_id,
                                                                     protein_length,
                                                                     CDS_exon_number,
                                                                     self.make_string_from_nested_list_of_ints(protein_segments_for_string),
                                                                     self.make_string_from_nested_list_of_ints(protein_segment_phases_list))
                    #print output_string

                    if translation != CDS_translation:
                        print("\tWARNING!!! Translational discrepancy in %s! Ignoring..." % protein_id)

                        print(output_full_string)
                        ignored_protein_counter += 1
                        continue
                    out_full_fd.write(output_full_string)
                    out_simple_fd.write(output_simple_string)

        print("\n\n%i proteins were ignored due to translational discrepancy" % ignored_protein_counter)

    @staticmethod
    def record_generator(annotations_dict, sequence_dict, feature_types_list):

        def reccursive_subfeature_retrival(feature):
            #print "aaaaaaaaaaaaaaaaaaaa"
            for subfeature in feature.sub_features:
                #print record_id, feature.id, subfeature.id, subfeature.type, feature.sub_features
                if subfeature.type in feature_types_list:
                    # print subfeature
                    sequence = subfeature.extract(sequence_dict[record_id].seq)
                    #record = SeqRecord(sequence, id=subfeature.id)
                    # print(record)
                    yield SeqRecord(sequence, id=subfeature.id, description=subfeature.qualifiers["Name"][0] \
                        if "Name" in subfeature.qualifiers else "")
                else:
                    for record in reccursive_subfeature_retrival(subfeature):
                        yield record

        for record_id in annotations_dict:
            for feature in annotations_dict[record_id].features:
                if feature.type in feature_types_list:
                    sequence = feature.extract(sequence_dict[record_id].seq)
                    #record = SeqRecord(sequence, id=feature.id)
                    #print(record)
                    yield SeqRecord(sequence, id=feature.id, description=feature.qualifiers["Name"][0] \
                          if "Name" in feature.qualifiers else "")

                    #print feature_types_list
                #print feature.sub_features
                for record in reccursive_subfeature_retrival(feature):
                    yield record

    @staticmethod
    def find_cds_coordinates_in_transcript_by_pep(transcript_dict, protein_dict, correspondence_dict,
                                                  verbose=True, genetic_code_table=1):

        transcript_with_no_protein_id_list = IdList()
        multiple_protein_hits_list = IdList()
        no_protein_hits_list = IdList()
        cds_coordinates_dict = OrderedDict()

        for transcript_id in transcript_dict:
            if transcript_id not in correspondence_dict:
                transcript_with_no_protein_id_list.append(transcript_id)
                if verbose:
                    print("Protein not found for transcript %s" % transcript_id)
                continue

            peptide = protein_dict[correspondence_dict[transcript_id]].seq

            reg_exp = re.compile(str(peptide))

            pep_frame_list = []
            matches_list = []

            for frame in (0, 1, 2):
                #pep_frame_list.append(transcript_dict[transcript_id].seq[frame:].translate(table=genetic_code_table))
                match_iterator = reg_exp.finditer(str(transcript_dict[transcript_id].seq[frame:].translate(table=genetic_code_table)))
                matches_coords_list = []
                for match in match_iterator:
                    #1-based coordinates for output
                    matches_coords_list.append((match.start() * 3 + frame + 1, match.end() * 3 + frame))

                if matches_coords_list:
                    matches_list += matches_coords_list

            if len(matches_list) == 0:
                print("No hits for protein %s in transcript %s. Skipping..." % (correspondence_dict[transcript_id],
                                                                                transcript_id))
                no_protein_hits_list.append(transcript_id)
                continue

            if len(matches_list) > 1:
                print("Multiple hits for protein %s in transcript %s. Skipping..." % (correspondence_dict[transcript_id],
                                                                                      transcript_id))
                multiple_protein_hits_list.append(transcript_id)
                continue
            #print matches_list
            cds_coordinates_dict[transcript_id] = matches_list[0]

        return cds_coordinates_dict, \
               transcript_with_no_protein_id_list, \
               no_protein_hits_list, \
               multiple_protein_hits_list

    def find_cds_coordinates_in_transcript_by_pep_from_file(self, transcript_file, protein_file, correspondence_file,
                                                            output_prefix, parsing_mode="parse", verbose=True,
                                                            format="fasta", transcript_index_file=None,
                                                            protein_index_file=None, genetic_code_table=1):
        transcript_dict = self.parse_seq_file(transcript_file, mode=parsing_mode, format=format,
                                              index_file=transcript_index_file)
        protein_dict = self.parse_seq_file(protein_file, mode=parsing_mode, format=format,
                                           index_file=protein_index_file)
        correspondence_dict = SynDict(filename=correspondence_file)

        cds_coordinates_dict, transcript_with_no_protein_id_list, no_protein_hits_list, multiple_protein_hits_list = \
            self.find_cds_coordinates_in_transcript_by_pep(transcript_dict, protein_dict, correspondence_dict,
                                                           verbose=verbose, genetic_code_table=genetic_code_table)

        cds_coordinates_file = "%s.cds.coordinates" % output_prefix
        transcript_with_no_protein_id_file = "%s.no_pep.ids" % output_prefix
        no_protein_hits_file = "%s.no_pep_hits.ids" % output_prefix
        multiple_protein_hits_file = "%s.mult_pep.ids" % output_prefix

        transcript_with_no_protein_id_list.write(transcript_with_no_protein_id_file)
        no_protein_hits_list.write(no_protein_hits_file)
        multiple_protein_hits_list.write(multiple_protein_hits_file)

        with open(cds_coordinates_file, "w") as coord_fd:
            for transcript_id in cds_coordinates_dict:
                coord_fd.write("%s\t%i\t%i\n" % (transcript_id, cds_coordinates_dict[transcript_id][0],
                                                 cds_coordinates_dict[transcript_id][1]))

    @staticmethod
    def split_ids_from_len_file_by_len(len_file, output_prefix, len_column=1, id_column=0):
        len_reverse_dict = SynDict(filename=len_file, allow_repeats_of_key=True,
                                   key_index=len_column, value_index=id_column, separator="\t")
        #print len_reverse_dict
        for length in len_reverse_dict:
            len_reverse_dict[length] = IdList(len_reverse_dict[length])
            out_file = "%s.%s.ids" % (output_prefix, length)
            len_reverse_dict[length].write(out_file)
        #print len_reverse_dict
        return len_reverse_dict

    @staticmethod
    def get_region_position_dict(location, segment_length):
        """
        location = [start, end) python notation
        """

        position_dict = OrderedDict()

        start_segment = int(location[0]/segment_length)
        end_segment = int((location[1] - 1)/segment_length)

        if start_segment == end_segment:
            position_dict[start_segment] = (location[0] % segment_length, ((location[1] - 1) % segment_length) + 1)
        elif end_segment - start_segment == 1:
            position_dict[start_segment] = (location[0] % segment_length, segment_length)
            position_dict[end_segment] = (0, ((location[1] - 1) % segment_length) + 1)
        else:
            position_dict[start_segment] = (location[0] % segment_length, segment_length)
            for i in range(start_segment + 1, end_segment):
                position_dict[i] = (0, segment_length)
            position_dict[end_segment] = (0, ((location[1] - 1) % segment_length) + 1)

        return position_dict

    def draw_string_regions(self, sequence, location_list, symbol_list,
                            overlap_symbol="#", line_per_record=False, segment_length=120,
                            num_of_spaces=3, num_of_space_lines=1, empty_symbol=" "):
        """
        TODO: add check for overlaps
        """

        sequence_length = len(sequence)

        last_segment_len = sequence_length % segment_length
        number_of_full_segments = (sequence_length - last_segment_len) / segment_length

        location_dict_list = []

        for location in location_list:
            location_dict_list.append(self.get_region_position_dict(location, segment_length))

        output_string = ""

        #print number_of_full_segments

        for i in range(0, (number_of_full_segments + 1) if last_segment_len > 1 else number_of_full_segments):
            output_string += "%12i%s%s%s%12i\n" % (i * segment_length + 1,
                                                   " " * num_of_spaces,
                                                   sequence[i * segment_length: (i + 1) * segment_length],
                                                   " " * num_of_spaces,
                                                   sequence_length if i == number_of_full_segments else ((i + 1) * segment_length))

            #regions_list = []

            #TODO: below is very stupid realization of possible overlap case, redo using analytical calculations

            region_string = " " * last_segment_len if i == number_of_full_segments else " " * segment_length

            for location_dict, symbol in zip(location_dict_list, symbol_list):
                if i in location_dict:
                    replacement = ""
                    for pos_in_segment in range(location_dict[i][0], location_dict[i][1]):
                        replacement += symbol if region_string[pos_in_segment] == empty_symbol else overlap_symbol
                    region_string = region_string[0:location_dict[i][0]] + replacement + region_string[location_dict[i][1]:]

            output_string += "%s%s%s%s%s\n" % (" " * 12,
                                               " " * num_of_spaces,
                                               region_string,
                                               " " * num_of_spaces,
                                               " " * 12)

            output_string += "\n" * num_of_space_lines

        return output_string

    def prepare_reference_for_GATK(self, reference, picard_dir="", samtools_dir=""):
        from Tools.Samtools import SamtoolsV1
        from Tools.Picard import CreateSequenceDictionary

        SamtoolsV1.path = samtools_dir
        CreateSequenceDictionary.jar_path = picard_dir

        SamtoolsV1.check_for_fasta_index(reference)
        CreateSequenceDictionary.check_for_fasta_dict(reference)

    def count_softmasked_nucleotides(self, record_dict, verbose=False, stats_file=None):

        total_softmasked_nucleotides = 0
        total_length = 0
        softmasked_nucleotides_dict = SynDict()
        percent_softmasked_nucleotides_dict = SynDict()

        for record_id in record_dict:
            record_len = len(record_dict[record_id].seq)
            softmasked_nucleotides_dict[record_id] = sum(c.islower() for c in record_dict[record_id].seq)
            total_softmasked_nucleotides += softmasked_nucleotides_dict[record_id]
            total_length += record_len

            percent_softmasked_nucleotides_dict[record_id] = float(softmasked_nucleotides_dict[record_id]) / float(record_len)

        stats_string = "Total length: %i\n" % total_length
        stats_string += "Softmasked nucleotides: %i\n" % total_softmasked_nucleotides
        stats_string += "Softmasked: %.2f\n" % (float(total_softmasked_nucleotides) / float(total_length))

        if verbose:
            print(stats_string)

        if stats_file:
            with self.metaopen(stats_file, "w") as out_fd:
                out_fd.write(stats_string)

        return softmasked_nucleotides_dict, percent_softmasked_nucleotides_dict

    def count_softmasked_nucleotides_from_file(self, sequence_file, output_prefix, verbose=False, parsing_mode="parse",
                                               format="fasta", index_file=None):

        record_dict = self.parse_seq_file(sequence_file, parsing_mode, format, index_file=index_file)

        softmasked_counts_file = "%s.sofmasking.counts" % output_prefix
        softmasked_percent_file = "%s.sofmasking.persent" % output_prefix
        softmasked_stats_file = "%s.sofmasking.stats" % output_prefix

        softmasked_nucleotides_dict, percent_softmasked_nucleotides_dict = self.count_softmasked_nucleotides(record_dict,
                                                                                                             verbose=verbose,
                                                                                                             stats_file=softmasked_stats_file)
        softmasked_nucleotides_dict.write(softmasked_counts_file)
        percent_softmasked_nucleotides_dict.write(softmasked_percent_file)

        return softmasked_nucleotides_dict, percent_softmasked_nucleotides_dict

    @staticmethod
    def generated_splited_regions(regions_list, sequence_len, min_length=1):
        sorted_region_list = sorted(regions_list)

        if sorted_region_list[-1][-1] > sequence_len:
            raise ValueError("ERROR!!! Last region is out of sequence!")

        splited_regions = []

        prev_end = 0

        for region in regions_list:
            if region[0] - prev_end >= min_length:
                splited_regions.append([prev_end, region[0]])
            prev_end = region[1]

        if sequence_len - prev_end >= min_length:
            splited_regions.append([prev_end, sequence_len])

        return splited_regions

    def split_sequence_by_regions(self, sequence_dict, regions_dict, retain_description=False,
                                  min_length=1, output_prefix=None):

        splited_sequence_dict = OrderedDict()

        splited_seq_ids = IdList()
        skipped_seq_ids = IdList()
        unchanged_seq_ids = IdList()

        retained_segments_dict = OrderedDict()

        for seq_id in sequence_dict:
            if seq_id not in regions_dict:
                splited_sequence_dict[seq_id] = sequence_dict[seq_id]
                unchanged_seq_ids.append(seq_id)
                continue

            splited_regions = self.generated_splited_regions(regions_dict[seq_id],
                                                             sequence_len=len(sequence_dict[seq_id].seq),
                                                             min_length=min_length)
            if not splited_regions:
                skipped_seq_ids.append(seq_id)
                continue

            splited_seq_ids.append(seq_id)
            retained_segments_dict[seq_id] = deepcopy(splited_regions)

            for i in range(0, len(splited_regions)):

                region_id = "%s_%i" % (seq_id, i + 1)
                splited_sequence_dict[region_id] = SeqRecord(seq=sequence_dict[seq_id].seq[splited_regions[i][0]:splited_regions[i][1]],
                                                             id=region_id,
                                                             description=sequence_dict[seq_id].description if retain_description else "")

        if output_prefix:
            splited_seq_ids.write("%s.splited_seqs.ids" % output_prefix)
            skipped_seq_ids.write("%s.skipped_seqs.ids" % output_prefix)
            unchanged_seq_ids.write("%s.unchanged_seqs.ids" % output_prefix)

            with self.metaopen("%s.retained_segments.bed" % output_prefix, "w") as out_fd:
                for seq_id in retained_segments_dict:
                    for region in retained_segments_dict[seq_id]:
                        out_fd.write("\t".join(map(str,
                                                   [seq_id,
                                                    region[0] + 1,
                                                    region[1]])) + "\n")


        return splited_sequence_dict

    def split_sequence_by_regions_from_file(self, sequence_file, regions_file, output_prefix,
                                            retain_description=False,
                                            min_length=1, parsing_mode="parse",
                                            scaffold_column_index=0,
                                            start_column_index=1,
                                            end_column_index=2,
                                            coordinates_type="1-based",
                                            input_separator="\t",
                                            sequence_format="fasta"):

        from Routines import AnnotationsRoutines

        sequence_dict = self.parse_seq_file(sequence_file, parsing_mode, format=sequence_format)

        regions_dict = AnnotationsRoutines.merge_overlapping_feature_in_simple_format(regions_file,
                                                                                      scaffold_column_index,
                                                                                      start_column_index,
                                                                                      end_column_index,
                                                                                      output_file=None,
                                                                                      output_separator="\t",
                                                                                      comments_prefix="#",
                                                                                      input_separator=input_separator,
                                                                                      coordinates_type=coordinates_type,
                                                                                      return_seqfeature_dict=False,
                                                                                      feature_type=None)

        splited_sequence_dict = self.split_sequence_by_regions(sequence_dict, regions_dict,
                                                               retain_description=retain_description,
                                                               min_length=min_length,
                                                               output_prefix=output_prefix,
                                                               )

        SeqIO.write(self.record_by_expression_generator(splited_sequence_dict, id_file=None),
                    "%s%s" % (output_prefix, self.split_filename(sequence_file)[-1]),
                    format=sequence_format)


def get_lengths(record_dict, out_file="lengths.t", write=False, write_header=True):
    lengths_dict = SynDict()
    lengths_dict.header = "#record\tlength\n"
    for record_id in record_dict:
        lengths_dict[record_id] = len(record_dict[record_id])

    if write:
        lengths_dict.write(out_file, header=write_header)

    return lengths_dict


def get_feature_lengths(record_dict):
    lengths_dict = OrderedDict({})

    for record_id in record_dict:
        lengths_dict[record_id] = {}
        for feature in record_dict[record_id].features:
            if feature.type not in lengths_dict[record_id]:
                lengths_dict[record_id][feature.type] = []
            lengths_dict[record_id][feature.type].append(len(feature))
            if feature.type == "gene":
                tmp_dict = {}
                for sub_feature in feature.sub_features:
                    if sub_feature.type not in tmp_dict:
                        tmp_dict[sub_feature.type] = 0
                    if sub_feature.type not in lengths_dict[record_id]:
                        lengths_dict[record_id][sub_feature.type] = []
                    tmp_dict[sub_feature.type] += len(sub_feature)

                for sub_feature_type in tmp_dict:
                    lengths_dict[record_id][sub_feature_type].append(tmp_dict[sub_feature_type])
    record_ids = lengths_dict.keys()
    lengths_dict["all"] = {}
    for record_id in record_ids:
        for feature_type in lengths_dict[record_id]:
            if feature_type not in lengths_dict["all"]:
                lengths_dict["all"][feature_type] = []
            lengths_dict["all"][feature_type] += lengths_dict[record_id][feature_type]

    #for record_id in lengths_dict:
    #    for feature_type in lengths_dict[record_id]:
    #        lengths_dict[record_id][feature_type] = np.array(lengths_dict[record_id][feature_type], dtype=int)
    return lengths_dict


def feature_lengths_collapse_records(lengths_dict,
                                     synonym_dict=None):
    tmp_dict = synonym_dict if synonym_dict is not None else {}
    collapsed_dict = {}
    for record_id in lengths_dict:
        if record_id == "all":
            continue
        for feature_type in lengths_dict[record_id]:
            tmp_type = tmp_dict[feature_type] if feature_type in tmp_dict else feature_type
            if tmp_type not in collapsed_dict:
                collapsed_dict[tmp_type] = lengths_dict[record_id][feature_type]
            else:
                collapsed_dict[tmp_type] += lengths_dict[record_id][feature_type]

    return collapsed_dict


def get_total_feature_lengths(lengths_dict, out_filename=None):
    total_lengths = TwoLvlDict({})
    for record_id in lengths_dict:
        total_lengths[record_id] = {}
        for feature_type in lengths_dict[record_id]:
            total_lengths[record_id][feature_type] = sum(lengths_dict[record_id][feature_type])
    if out_filename is not None:
        total_lengths.write(out_filename, absent_symbol="0")

    return total_lengths


# --------------------Search--------------------------
#              !!!!IMPORTANT!!!!!!
# in this section python notation for coordinates inside sequence is used
def find_homopolymer_end(seq, nucleotide, seq_length, start, search_type="perfect",
                       max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2):
    shift_list = []
    number_of_insertions = 0
    total_insertion_length = 0
    insertion_length = 0
    i = start + 1
    if search_type == 'perfect':
        while i < seq_length:
            if seq[i] != nucleotide:
                return i, None
            i += 1
        return seq_length, None
    else:
        while i < seq_length:
            if seq[i] != nucleotide:
                if seq[i-1] == nucleotide:
                    shift_list.append(i)
                    insertion_length = 1
                else:
                    insertion_length += 1
                number_of_insertions += 1
                total_insertion_length += 1
                if number_of_insertions > max_number_of_insertions or insertion_length > max_single_insert_size:
                    end = shift_list[-1]
                    break
                if max_total_insert_length:
                    if total_insertion_length > max_total_insert_length:
                        end = shift_list[-1]
                        break
            i += 1
        else:
            if seq[-1] == nucleotide:
                end = seq_length
            else:
                end = shift_list[0]

        new_start = shift_list[0] if shift_list else end
        return end, new_start




# ----------------------Filters-----------------------


def filter_sequences(record_dict, expression):
    record_id_list = []
    for record_id in record_dict:
        if expression(record_dict[record_id]):
            record_id_list.append(record_id)
    return record_id_list

# --------------------Generators----------------------


def rev_com_generator(record_dict, yield_original_record=False, rc_suffix="_rc", add_rc_suffix="False"):
    for record_id in record_dict:
        reverse_seq = record_dict[record_id].seq.reverse_complement()
        record = deepcopy(record_dict[record_id])
        record.seq = reverse_seq
        record.id = record_id + rc_suffix if yield_original_record or add_rc_suffix else record.id
        if yield_original_record:
            yield record_dict[record_id]
        yield record


def find_gaps(record_dict):
    gap_reg_exp = re.compile("N+", re.IGNORECASE)
    gaps_dict = {}
    for region in record_dict:
        gaps_dict[region] = SeqRecord(seq=record_dict[region].seq,
                                      id=record_dict[region].id,
                                      description=record_dict[region].description)
        gaps = gap_reg_exp.finditer(str(record_dict[region].seq))  # iterator with
        for match in gaps:
            gaps_dict[region].features.append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                              type="gap", strand=None))
    return gaps_dict

# legacy function, moved to class SequenceRoutines as method, possibly all usages were changed
"""

def record_by_id_generator(record_dict, id_list):
    for record_id in id_list:
        if record_id in record_dict:
            yield record_dict[record_id]
        else:
            sys.stderr.write("Not found: %s\n" % record_id)
"""


def record_by_expression_generator(record_dict, expression=None, id_file="passed_records.ids"):
    """
    :param record_dict: dictionary containing Biopython SeqRecords as values
    :param expression: function to apply to all records in record_dict. If it returns True record will be yielded
    :param id_file: file to write ids of records passed expression
    :return: None
    """
    id_fd = open(id_file, "w")
    if expression is None:
        for record_id in record_dict:
            id_fd.write(record_id + "\n")
            yield record_dict[record_id]
    else:
        for record_id in record_dict:
            if expression(record_dict[record_id]):
                id_fd.write(record_id + "\n")
                yield record_dict[record_id]
    id_fd.close()




def record_generator(annotations_dict, sequence_dict, feature_types_list):
        for record_id in annotations_dict:
            for feature in annotations_dict[record_id].features:
                if feature.type in feature_types_list:
                    sequence = feature.extract(sequence_dict[record_id].seq)
                    #record = SeqRecord(sequence, id=feature.id)
                    #print(record)
                    yield SeqRecord(sequence, id=feature.id, description=feature.qualifiers["Name"][0] \
                          if "Name" in feature.qualifiers else "")

def get_kmer_dict(sequence, kmer_length, start=1, end=None):  # positions are one-based
    kmer_dict = OrderedDict()
    length = len(sequence)
    if end > length:
        raise ValueError("End position should be less than length")
    stop = length if end is None else end
    for i in range(start-1, stop-kmer_length):
        kmer_dict[i+1] = sequence[i:i+kmer_length]
    return kmer_dict


def get_kmer_dict_as_seq_records(sequence, kmer_length, start=1, end=None, id_prefix="id"):  # positions are one-based
    kmer_dict = OrderedDict()
    length = len(sequence)
    if end > length:
        raise ValueError("End position should be less than length")
    stop = length if end is None else end
    for i in range(start-1, stop-kmer_length):
        record_id = "%s_%i-%i" % (id_prefix, i + 1, i + kmer_length)
        kmer_dict[record_id] = SeqRecord(seq=sequence[i:i+kmer_length], id=record_id)
    return kmer_dict


def split_record_ids_by_expression(sequence_dict, expression):
    """
    splits record ids based on parameter of record
    """
    id_dict = OrderedDict()
    for record_id in sequence_dict:

        value = expression(sequence_dict[record_id])
        if value not in id_dict:
            id_dict[value] = [record_id]
        else:
            id_dict[value].append(record_id)

    return id_dict
