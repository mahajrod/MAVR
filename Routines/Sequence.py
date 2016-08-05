__author__ = 'mahajrod'
import os
import re
import sys
import pickle

from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from CustomCollections.GeneralCollections import TwoLvlDict, SynDict, IdList, IdSet
from Routines.Functions import output_dict


class SequenceRoutines():

    def __init__(self):
        pass

    @staticmethod
    def get_lengths(record_dict, out_file=None, close_after_if_file_object=False):
        lengths_dict = SynDict()
        for record_id in record_dict:
            lengths_dict[record_id] = len(record_dict[record_id])

        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)

        return lengths_dict

    @staticmethod
    def get_lengths_from_seq_file(input_file_list, format="fasta", out_file=None, close_after_if_file_object=False):
        record_dict = SeqIO.index_db("tmp.idx", input_file_list, format=format)
        lengths_dict = SynDict()
        for record_id in record_dict:
            lengths_dict[record_id] = len(record_dict[record_id])

        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)

        os.remove("tmp.idx")
        return lengths_dict

    @staticmethod
    def record_by_id_generator(record_dict, id_list, verbose=False):
        for record_id in id_list:
            if record_id in record_dict:
                yield record_dict[record_id]
            else:
                if verbose:
                    sys.stderr.write("Not found: %s\n" % record_id)
                    
    @staticmethod
    def record_from_dict_generator(record_dict):
        for record_id in record_dict:
            yield record_dict[record_id]

    def extract_sequence_by_ids(self, sequence_file, id_file, output_file, format="fasta", verbose=False):
        tmp_index_file = "tmp.idx"
        id_list = IdList()
        id_list.read(id_file)
        if verbose:
            print("Parsing %s..." % (sequence_file if isinstance(id_file, str) else ",".join(id_file)))

        sequence_dict = SeqIO.index_db(tmp_index_file, sequence_file, format=format)
        SeqIO.write(self.record_by_id_generator(sequence_dict, id_list, verbose=verbose),
                    output_file, format=format)
        os.remove(tmp_index_file)

    def extract_sequences_by_length_from_file(self, input_file, output_file, min_len=1, max_len=None, format="fasta",
                                              tmp_index_file="tmp.idx", id_file=None):

        if (min_len is None) and (max_len is None):
            raise ValueError("Both minimum and maximum lengths were not set")
        elif (min_len is not None) and (max_len is not None) and (min_len > max_len):
            raise ValueError("Minimum length is greater then maximum lengths")
        elif (min_len is not None) and (min_len < 0):
            raise ValueError("Minimum length is below zero")
        elif (max_len is not None) and (max_len < 0):
            raise ValueError("Maximum length is below zero")

        sequence_dict = SeqIO.index_db(tmp_index_file, input_file, format=format)

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
                            print "BBBBBBBBBB"
                            if expression(record_dict[record_id]):
                                print("CCCCCCCC")
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
    def parse_seq_file(input_file, mode, format="fasta", index_file=None):
        if mode == "index_db":
            index = index_file if index_file else "tmp.idx"
            seq_dict = SeqIO.index_db(index, input_file, format=format)
        elif mode == "index":
            seq_dict = SeqIO.index(input_file, format=format)
        elif mode == "parse":
            seq_dict = SeqIO.to_dict(SeqIO.parse(input_file, format=format))

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
                                      prefix_of_file_inframe_stop_codons_seqs="cds_with_in_frame_stop_codons"):
        cds_dict = SeqIO.index_db("tmp.idx", input_file, format=format)
        SeqIO.write(self.translated_seq_generator(cds_dict, id_expression=id_expression,
                                                  genetic_code_table=genetic_code_table,
                                                  translate_to_stop=translate_to_stop,
                                                  prefix_of_file_inframe_stop_codons_seqs=prefix_of_file_inframe_stop_codons_seqs),
                    output_file, format=format)
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

    def check_proteins_for_stop_codons_from_file(self, pep_file, output_prefix, stop_codon_symbol_set=("*", "."), format="fasta"):
        pep_with_stop_codons_ids_file = "%s.ids" % output_prefix
        pep_with_stop_codons_seq_file = "%s.pep" % output_prefix
        pep_dict = SeqIO.index_db("tmp.idx", pep_file, format=format)

        pep_with_stop_codons_ids = self.check_proteins_for_stop_codons(pep_dict,
                                                                       stop_codon_symbol_set=stop_codon_symbol_set)
        pep_with_stop_codons_ids.write(pep_with_stop_codons_ids_file)
        SeqIO.write(self.record_by_id_generator(pep_dict, pep_with_stop_codons_ids), pep_with_stop_codons_seq_file,
                    format=format)
        os.remove("tmp.idx")
        return pep_with_stop_codons_ids

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

        for cds_id in cds_dict:
            cds_pep = cds_dict[cds_id].seq.translate(to_stop=True, table=genetic_code_table)
            for pep_id in pep_dict:
                if include_id_check:
                    if cds_id in pep_dict:
                        cds_pep_accordance_dict[cds_id] = cds_id
                        if parsing_mode == "parse":
                            pep_dict.pop(pep_id, None)
                            break

                if cds_pep == pep_dict[pep_id].seq:
                    cds_pep_accordance_dict[cds_id] = pep_id
                    if parsing_mode == "parse":
                        pep_dict.pop(pep_id, None)
                    break
            else:
                if verbose:
                    print("Protein was not found for %s CDS" % cds_id)

        return cds_pep_accordance_dict

    def get_cds_to_pep_accordance_from_files(self, cds_file, pep_file, output_file, format="fasta",
                                             verbose=True, parsing_mode="parse", index_file_suffix="tmp.idx",
                                             genetic_code_table=1, include_id_check=False):

        cds_dict = self.parse_seq_file(cds_file, mode=parsing_mode, format=format,
                                       index_file="cds_%s" % index_file_suffix)
        pep_dict = self.parse_seq_file(pep_file, mode=parsing_mode, format=format,
                                       index_file="pep_%s" % index_file_suffix)

        cds_pep_accordance_dict = self.get_cds_to_pep_accordance(cds_dict, pep_dict, verbose=verbose,
                                                                 parsing_mode=parsing_mode,
                                                                 genetic_code_table=genetic_code_table,
                                                                 include_id_check=include_id_check)
        cds_pep_accordance_dict.write(output_file)

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
    def renamed_records_generator(record_dict, syn_dict=None, expression=None, clear_description=False):
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

            yield record

    def rename_records_from_files(self, input_file, output_file, synonyms_file=None, format="fasta", header=False,
                                  separator="\t", key_index=0, value_index=1, syn_expression=None, comments_prefix=None,
                                  clear_description=False, record_id_expression=None):
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

        record_dict = SeqIO.index_db("temp.idx", input_file, format=format)

        SeqIO.write(self.renamed_records_generator(record_dict, syn_dict=syn_dict, expression=record_id_expression,
                                                   clear_description=clear_description),
                    output_file, format=format)

        os.remove("temp.idx")

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
    def make_region_bed_file(record_dict, output_file, white_list=None, black_list=None, output_format="0-based"):
        with open(output_file, "w") as out_fd:
            for record_id in record_dict:
                if white_list:
                    if record_id not in white_list:
                        continue
                if black_list:
                    if record_id in black_list:
                        continue
                bed_string = "%s" % record_id
                if output_format == "0-based":
                    bed_string += "\t0\t%i\n" % len(record_dict[record_id])
                out_fd.write(bed_string)

    def make_region_bed_file_from_file(self, seq_file, output_file, white_id_file=None, black_id_file=None,
                                       output_format="0-based", input_format="fasta"):
        record_dict = SeqIO.index_db("tmp.idx", seq_file, format=input_format)
        black_id_list = IdList()
        white_id_list = IdList()
        if black_id_file:
            black_id_list.read(black_id_file)
        if white_id_file:
            white_id_list.read(white_id_file)

        self.make_region_bed_file(record_dict, output_file, white_list=white_id_list, black_list=black_id_list,
                                  output_format=output_format)

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


def get_lengths(record_dict, out_file="lengths.t", write=False, write_header=True):
    lengths_dict = OrderedDict({})
    for record_id in record_dict:
        lengths_dict[record_id] = len(record_dict[record_id])

    if write:
        output_dict(lengths_dict, out_file=out_file, write=write,
                    header_tuple=("record", "length") if write_header else None)

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
