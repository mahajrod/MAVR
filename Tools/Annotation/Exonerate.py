#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from multiprocessing import Pool
from collections import OrderedDict
from Tools.Abstract import Tool

from CustomCollections.GeneralCollections import SynDict, IdList


from Bio import SeqIO

class VulgarAlignment:
    def __init__(self, exonerate_vulgar_string):
        string_list = exonerate_vulgar_string.split()
        self.query_id = string_list[0]
        self.query_start = string_list[1]
        self.query_end = string_list[2]
        self.query_strand = string_list[3]
        self.target_id = string_list[4]
        self.target_start = string_list[5]
        self.target_end = string_list[6]
        self.target_strand = string_list[7]
        self.alignment_score = string_list[8]
        self.vulgar_dict = self.parse_vulgar_string(string_list[9:])

    @staticmethod
    def parse_vulgar_string(vulgar_string):
        vulgar_dict = OrderedDict()
        if isinstance(vulgar_string, str):
            vulgar_list = vulgar_string.strip().split()
        else:
            vulgar_list = vulgar_string
        if len(vulgar_list) % 3 != 0:
            raise ValueError("ERROE: Input vulgar string is malformed!")
        for i in range(0, int(len(vulgar_list)/3)):
            vulgar_dict[vulgar_list[i*3]] = (vulgar_list[i*3+1], vulgar_list[i*3+2])
        return vulgar_list

    def parse_exonerate_vulgar_string(self, exonerate_vulgar_string):
        string_list = exonerate_vulgar_string.split()
        query_id = string_list[0]
        query_start = string_list[1]
        query_end = string_list[2]
        query_strand = string_list[3]
        target_id = string_list[4]
        target_start = string_list[5]
        target_end = string_list[6]
        target_strand = string_list[7]
        alignment_score = string_list[8]
        vulgar_dict = self.parse_vulgar_string(string_list[9:])

        return ((query_id, query_start, query_end, query_strand), (target_id, target_start, target_end, target_strand),
                alignment_score, vulgar_dict)


class Exonerate(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "exonerate", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(model, query_type=None, target_type=None,
                             show_alignment=None, show_sugar=True, show_cigar=None,
                             show_vulgar=None, show_query_gff=None, show_target_gff=None,
                             number_of_results_to_report=None, other_options=None, annotation_file=None,
                             softmasked_target=False, softmasked_query=False):

        options = " --model %s" % model
        options += " --showalignment" if show_alignment else ""
        options += " --showsugar" if show_sugar else ""
        options += " --showcigar" if show_cigar else ""
        options += " --showvulgar" if show_vulgar else ""
        options += " --showquerygff" if show_query_gff else ""
        options += " --showtargetgff" if show_target_gff else ""
        options += " --softmasktarget" if softmasked_target else ""
        options += " --softmaskquery" if softmasked_query else ""
        options += " -Q %s" % query_type if query_type else ""
        options += " -T %s" % target_type if target_type else ""
        options += " -n %i" % number_of_results_to_report if number_of_results_to_report else ""
        options += " %s" % other_options if other_options else ""
        options += " --annotation %s" % annotation_file if annotation_file else ""

        return options

    def parallel_alignment(self, query_file, target_file, model, num_of_recs_per_file=None,
                           show_alignment=None, show_sugar=True, show_cigar=None,
                           show_vulgar=None, show_query_gff=None, show_target_gff=None,
                           store_intermediate_files=True,
                           annotation_file=None,
                           splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                           number_of_results_to_report=None,
                           other_options=None,
                           num_of_files=None,
                           converted_output_dir="converted_output", parsing_mode="parse", index_file=None,
                           external_process_pool=None,
                           softmasked_target=False, softmasked_query=False):
        splited_filename = self.split_filename(query_file)
        self.split_fasta(query_file, splited_fasta_dir, num_of_recs_per_file=num_of_recs_per_file,
                         num_of_files=num_of_files,
                         output_prefix=splited_filename[1],
                         parsing_mode=parsing_mode, index_file=index_file)

        common_options = self.parse_common_options(model, show_alignment=show_alignment,
                                                   show_sugar=show_sugar, show_cigar=show_cigar,
                                                   show_vulgar=show_vulgar, show_query_gff=show_query_gff,
                                                   show_target_gff=show_target_gff,
                                                   number_of_results_to_report=number_of_results_to_report,
                                                   other_options=other_options,
                                                   annotation_file=annotation_file,
                                                   softmasked_target=softmasked_target,
                                                   softmasked_query=softmasked_query)

        options_list = []
        splited_files = os.listdir(splited_fasta_dir)

        self.safe_mkdir(splited_result_dir)
        #save_mkdir(converted_output_dir)

        for filename in splited_files:
            filename_list = self.split_filename(filename)
            options = common_options
            options += " -q %s/%s" % (splited_fasta_dir, filename)
            options += " -t %s" % target_file
            options += " > %s/%s.output" % (splited_result_dir, filename_list[1])
            options_list.append(options)

        self.parallel_execute(options_list, external_process_pool=external_process_pool)

        if not store_intermediate_files:
            shutil.rmtree(splited_fasta_dir)
            #shutil.rmtree(splited_result_dir)
            #shutil.rmtree(converted_output_dir)

    def prepare_index(self, list_of_files, output_prefix, translated_index=None, memory_limit=1024):
        pass

    def parallel_alignment_server_mode(self, query_file, target_file, model, num_of_recs_per_file=None,
                                       show_alignment=None, show_sugar=True, show_cigar=None,
                                       show_vulgar=None, show_query_gff=None, show_target_gff=None,
                                       store_intermediate_files=True,
                                       annotation_file=None,
                                       splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                                       number_of_results_to_report=None,
                                       other_options=None,
                                       num_of_files=None,
                                       converted_output_dir="converted_output", parsing_mode="parse", index_file=None,
                                       external_process_pool=None):
        


        pass

    def prepare_data_for_target_alignment(self, query_fasta, target_fasta, correspondence_file, out_dir,
                                          correspondence_query_column=0, correspondence_target_column=1):

        query_dict = self.parse_seq_file(query_fasta, "parse")
        target_dict = self.parse_seq_file(target_fasta, "parse")

        self.safe_mkdir(out_dir)

        correspondence_dict = SynDict(filename=correspondence_file, allow_repeats_of_key=True,
                                      key_index=correspondence_query_column, value_index=correspondence_target_column)

        for query_id in correspondence_dict:
            query_outfile = "%s/%s.query.fasta" % (out_dir, query_id)
            target_outfile = "%s/%s.target.fasta" % (out_dir, query_id)

            SeqIO.write(self.record_by_id_generator(query_dict, [query_id]), query_outfile, format="fasta")
            SeqIO.write(self.record_by_id_generator(target_dict, correspondence_dict[query_id]),
                        target_outfile, format="fasta")

        queries_with_targets_set = set(correspondence_dict.keys())
        queries_set = set(query_dict.keys())

        return queries_with_targets_set, queries_set - queries_with_targets_set
    """
    def parallel_target_alignment(self, query_file, target_file, model, num_of_recs_per_file=None,
                                   show_alignment=None, show_sugar=True, show_cigar=None,
                                   show_vulgar=None, show_query_gff=None, show_target_gff=None,
                                   store_intermediate_files=True,
                                   annotation_file=None,
                                   splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                                   number_of_results_to_report=None,
                                   other_options=None,
                                   num_of_files=None,
                                   converted_output_dir="converted_output", parsing_mode="parse", index_file=None):
        splited_filename = self.split_filename(query_file)
        self.split_fasta(query_file, splited_fasta_dir, num_of_recs_per_file=num_of_recs_per_file,
                         num_of_files=num_of_files,
                         output_prefix=splited_filename[1],
                         parsing_mode=parsing_mode, index_file=index_file)

        common_options = self.parse_common_options(model, show_alignment=show_alignment,
                                                   show_sugar=show_sugar, show_cigar=show_cigar,
                                                   show_vulgar=show_vulgar, show_query_gff=show_query_gff,
                                                   show_target_gff=show_target_gff,
                                                   number_of_results_to_report=number_of_results_to_report,
                                                   other_options=other_options,
                                                   annotation_file=annotation_file)

        options_list = []
        splited_files = os.listdir(splited_fasta_dir)

        self.safe_mkdir(splited_result_dir)
        #save_mkdir(converted_output_dir)

        for filename in splited_files:
            filename_list = self.split_filename(filename)
            options = common_options
            options += " -q %s/%s" % (splited_fasta_dir, filename)
            options += " -t %s" % target_file
            options += " > %s/%s.output" % (splited_result_dir, filename_list[1])
            options_list.append(options)

        self.parallel_execute(options_list)

        if not store_intermediate_files:
            shutil.rmtree(splited_fasta_dir)
            #shutil.rmtree(splited_result_dir)
            #shutil.rmtree(converted_output_dir)

    """

    def split_output(self, exonerate_output_files, output_prefix, reference_protein_file, gene_prefix="GEN", transcript_prefix="TR", number_len=8,
                     ):
        """
        full_length hits - whole query was aligned
        """
        
        names_dict = {
                      "vulgar": "%s.vulgar" % output_prefix,
                      "cigar": "%s.cigar" % output_prefix,
                      "sugar": "%s.sugar" % output_prefix,
                      "alignment": "%s.alignment" % output_prefix,
                      "gff": "%s.gff" % output_prefix,
                      "target_gff": "%s.target.gff" % output_prefix,
                      "query_gff": "%s.query.gff" % output_prefix,

                      "splice": "%s.splice.gff" % output_prefix,
                      "exon": "%s.exon.gff" % output_prefix,
                      "intron": "%s.intron.gff" % output_prefix,
                      "cds": "%s.cds.gff" % output_prefix,
                      "gene": "%s.gene.gff" % output_prefix,
            
                      "vulgar_full_length_top": "%s.full_length_top.vulgar" % output_prefix,
                      "cigar_full_length_top": "%s.full_length_top.cigar" % output_prefix,
                      "sugar_full_length_top": "%s.full_length_top.sugar" % output_prefix,
                      "alignment_full_length_top": "%s.full_length_top.alignment" % output_prefix,
                      "gff_full_length_top": "%s.full_length_top.gff" % output_prefix,
                      "target_gff_full_length_top": "%s.full_length_top.target.gff" % output_prefix,
                      "query_gff_full_length_top": "%s.full_length_top.query.gff" % output_prefix,
            
                      "vulgar_other_top": "%s.other_top.vulgar" % output_prefix,
                      "cigar_other_top": "%s.other_top.cigar" % output_prefix,
                      "sugar_other_top": "%s.other_top.sugar" % output_prefix,
                      "alignment_other_top": "%s.other_top.alignment" % output_prefix,
                      "gff_other_top": "%s.other_top.gff" % output_prefix,
                      "target_gff_other_top": "%s.other_top.target.gff" % output_prefix,
                      "query_gff_other_top": "%s.other_top.query.gff" % output_prefix,
            
                      "vulgar_full_length_secondary": "%s.full_length_secondary.vulgar" % output_prefix,
                      "cigar_full_length_secondary": "%s.full_length_secondary.cigar" % output_prefix,
                      "sugar_full_length_secondary": "%s.full_length_secondary.sugar" % output_prefix,
                      "alignment_full_length_secondary": "%s.full_length_secondary.alignment" % output_prefix,
                      "gff_full_length_secondary": "%s.full_length_secondary.gff" % output_prefix,
                      "target_gff_full_length_secondary": "%s.full_length_secondary.target.gff" % output_prefix,
                      "query_gff_full_length_secondary": "%s.full_length_secondary.query.gff" % output_prefix,
            
                      "vulgar_other_secondary": "%s.other_secondary.vulgar" % output_prefix,
                      "cigar_other_secondary": "%s.other_secondary.cigar" % output_prefix,
                      "sugar_other_secondary": "%s.other_secondary.sugar" % output_prefix,
                      "alignment_other_secondary": "%s.other_secondary.alignment" % output_prefix,
                      "gff_other_secondary": "%s.other_secondary.gff" % output_prefix,
                      "target_gff_other_secondary": "%s.other_secondary.target.gff" % output_prefix,
                      "query_gff_other_secondary": "%s.other_secondary.query.gff" % output_prefix,

                      "stats_full_length_top": "%s.full_length_top.stats" % output_prefix,
                      "stats_other_top": "%s.other_top.stats" % output_prefix,
                      "stats_full_length_secondary": "%s.full_length_secondary.stats" % output_prefix,
                      "stats_other_secondary": "%s.other_secondary.stats" % output_prefix,
                      "stats": "%s.stats" % output_prefix,
                      "stats_hit": "%s.hit.stats" % output_prefix,
                      "stats_general": "%s.general.stats" % output_prefix,

                      "top_hits_gff" : "%s.top_hits.gff" % output_prefix,
                      "top_hits_vulgar" : "%s.top_hits.vulgar" % output_prefix,
                      "top_hits_sugar" : "%s.top_hits.sugar" % output_prefix,
                      "top_hits_target_gff" : "%s.top_hits.target.gff" % output_prefix,
                      "top_hits_query_gff" : "%s.top_hits.query.gff" % output_prefix,
                      "top_hits_simple" : "%s.top_hits.query.simple" % output_prefix,
                    }
        
        reference_protein_dict = self.parse_seq_file(reference_protein_file, "parse", format="fasta")

        fd_dict = {}
        for output_type in names_dict:
            fd_dict[output_type] = open(names_dict[output_type], "w")
            
        for output_type_entry in "stats_full_length_top", "stats_other_top", \
                                 "stats_full_length_secondary", "stats_other_secondary", "stats":
            fd_dict[output_type_entry].write("#query_id\tquery_len\tscore\thit_len\tquery_start\tquery_end\tgene_id\n")

        fd_dict["stats_hit"].write("#query_id\ttotal_hits\tfull_length_top_hits\tother_top_hit\tfull_length_secondary_hits\tother_secondary_hits\n")

        current_gene_index = 0
        gene_prefix = "%s%%0%ii" % (gene_prefix, number_len)
        transcript_prefix = "%s%%0%ii" % (transcript_prefix, number_len)
        full_length_top_counter = 0
        other_top_counter = 0
        full_length_secondary_counter = 0
        other_secondary_counter = 0
        previous_query_id = ""
        hit_counter = 0
        full_length_flag = False
        current_query_len = 0
        current_raw_score = 0
        current_query_start = 0
        current_query_end = 0
        current_hit_length = 0

        per_pep_full_length_top_hit_number = 0
        per_pep_other_top_hit_number = 0
        per_pep_full_length_secondary_hit_number = 0
        per_pep_other_secondary_hit_number = 0

        for filename in exonerate_output_files:
            index = 0

            with open(filename, "r") as in_fd:
                #print u
                #tmp = None

                for line in in_fd:
                    tmp = line

                    if tmp[:13] == "C4 Alignment:":
                        alignment_buffer = tmp

                        while True:
                            tmp = in_fd.next()
                            alignment_buffer += tmp

                            if "Query: " in tmp:
                                current_query_id = tmp.split()[1]
                                if current_query_id != previous_query_id:
                                    if previous_query_id != "":
                                        fd_dict["stats_hit"].write("%s\t%i\t%i\t%i\t%i\t%i\n" % (previous_query_id,
                                                                                                 hit_counter,
                                                                                                 per_pep_full_length_top_hit_number,
                                                                                                 per_pep_other_top_hit_number,
                                                                                                 per_pep_full_length_secondary_hit_number,
                                                                                                 per_pep_other_secondary_hit_number))

                                    previous_query_id = current_query_id
                                    hit_counter = 1

                                    per_pep_full_length_top_hit_number = 0
                                    per_pep_other_top_hit_number = 0
                                    per_pep_full_length_secondary_hit_number = 0
                                    per_pep_other_secondary_hit_number = 0

                                    current_query_len = len(reference_protein_dict[current_query_id].seq)
                                else:
                                    hit_counter += 1
                            elif "Raw score:" in tmp:
                                current_raw_score = tmp.strip().split()[-1]
                                #print current_raw_score
                            elif "Query range:" in tmp:
                                #print tmp
                                #print tmp.strip("Query range: ").split()[-1].split(" -> ")
                                current_query_start, current_query_end = map(int, tmp.strip().split("Query range: ")[-1].split(" -> "))
                                current_hit_length = current_query_end - current_query_start
                                
                                full_length_flag = True if current_hit_length == current_query_len else False
                                
                            elif "Target range:" in tmp:

                                if hit_counter == 1:
                                    if full_length_flag:
                                        output_type = "full_length_top"
                                        full_length_top_counter += 1
                                        per_pep_full_length_top_hit_number += 1
                                    else:
                                        other_top_counter += 1
                                        output_type = "other_top"
                                        per_pep_other_top_hit_number += 1
                                else:
                                    if full_length_flag:
                                        output_type = "full_length_secondary"
                                        full_length_secondary_counter += 1
                                        per_pep_full_length_secondary_hit_number += 1
                                    else:
                                        output_type = "other_secondary"
                                        other_secondary_counter += 1
                                        per_pep_other_secondary_hit_number += 1

                                fd_dict["alignment"].write(alignment_buffer)
                                fd_dict["alignment_" + output_type].write(alignment_buffer)
                                
                                break

                        while True:
                            tmp = next(in_fd, "")
                            if (tmp[0] == "c") or (tmp[0] == "v") or (tmp[0] == "s") or (tmp[0] == "#"):
                                break

                            fd_dict["alignment"].write(tmp)
                            fd_dict["alignment_" + output_type].write(tmp)
                            if tmp == "":
                                break
                    if tmp == "# --- START OF GFF DUMP ---\n":
                        fd_dict["gff"].write(tmp)
                        fd_dict["gff_" + output_type].write(tmp)
                        if index == 0:
                            fd_dict["query_gff"].write(tmp)
                            fd_dict["query_gff_" + output_type].write(tmp)
                        else:
                            fd_dict["target_gff"].write(tmp)
                            fd_dict["target_gff_" + output_type].write(tmp)
                        while True:
                            tmp = next(in_fd, "")
                            if ("\tgene\t" in tmp) and (index == 1):
                                #print tmp
                                current_gene_index += 1
                                #print current_gene_index, gene_prefix
                                current_gene_id = gene_prefix % current_gene_index
                                current_transcript_id = transcript_prefix % current_gene_index

                                hit_stat_str = "%s\t%s\t%s\t%i\t%i\t%i\t%s\n" % (current_query_id,
                                                                                 current_query_len,
                                                                                 current_raw_score,
                                                                                 current_hit_length,
                                                                                 current_query_start + 1,
                                                                                 current_query_end,
                                                                                 current_gene_id)
                                fd_dict["stats"].write(hit_stat_str)
                                fd_dict["stats_" + output_type].write(hit_stat_str)

                                line_list = tmp.strip().split("\t")
                                attr_list = map(lambda s: s.split(), line_list[-1].split(";"))
                                for i in range(0, len(attr_list)):
                                    if attr_list[i][0] == "gene_id":
                                        attr_list[i][1] = current_gene_id

                                line_list[-1] = ";".join(map(lambda s: " ".join(s), attr_list))
                                tmp = "\t".join(line_list) + "\n"

                                line_list[-1] = ("transcript_id %s ; " % current_transcript_id) + line_list[-1]
                                line_list[2] = "transcript"
                                transcript_line = "\t".join(line_list) + "\n"
                                #print current_gene_id
                                #print tmp
                                #print current_gene_id
                            elif (tmp[0] != "#") and (index == 1):
                                #print tmp
                                tmp_list = tmp.split("\t")

                                if ("\texon\t" in tmp) or ("\tcds\t" in tmp):
                                    tmp_list[-1] = ("gene_id %s ; transcript_id %s ; " % (current_gene_id, current_transcript_id)) + tmp_list[-1]
                                else:
                                    tmp_list[-1] = ("gene_id %s ; " % current_gene_id) + tmp_list[-1]
                                tmp = "\t".join(tmp_list)

                            fd_dict["gff"].write(tmp)
                            if index == 0:
                                fd_dict["query_gff"].write(tmp)
                                fd_dict["query_gff_" + output_type].write(tmp)
                            else:
                                fd_dict["target_gff"].write(tmp)
                                fd_dict["target_gff_" + output_type].write(tmp)
                                if "\tgene\t" in tmp:
                                    fd_dict["target_gff"].write(transcript_line)
                                    fd_dict["target_gff_" + output_type].write(transcript_line)
                            if tmp[0] != "#":
                                if "\tsplice" in tmp:
                                    fd_dict["splice"].write(tmp)
                                elif "\texon\t" in tmp:
                                    fd_dict["exon"].write(tmp)
                                elif "\tintron\t" in tmp:
                                    fd_dict["intron"].write(tmp)
                                elif "\tcds\t" in tmp:
                                    fd_dict["cds"].write(tmp)
                                elif "\tgene\t" in tmp:
                                    fd_dict["gene"].write(tmp)

                            if tmp == "# --- END OF GFF DUMP ---\n":
                                break
                        index = 1 if index == 0 else 0
                        continue

                    if tmp == "":
                        break
                    if tmp[:7] == "vulgar:":
                        fd_dict["vulgar"].write(tmp[7:])
                        fd_dict["vulgar_" + output_type].write(tmp[7:])
                    elif tmp[:6] == "sugar:":
                        fd_dict["sugar"].write(tmp[6:])
                        fd_dict["sugar_" + output_type].write(tmp[6:])
                    elif tmp[:6] == "cigar:":
                        fd_dict["cigar"].write(tmp[6:])
                        fd_dict["cigar" + output_type].write(tmp[6:])

        # extract top hits from vulgar
        awk_string_prefix = "awk '{if (substr($0,0,1) != \"#\") {if ($1 != SEQ_ID) {print $0}; SEQ_ID=$1}}' "
        os.system(awk_string_prefix + " %s > %s" % (names_dict["vulgar"], names_dict["top_hits_vulgar"]))
        os.system(awk_string_prefix + " %s > %s" % (names_dict["sugar"], names_dict["top_hits_sugar"]))
        os.system(awk_string_prefix + " %s > %s" % (names_dict["query_gff"], names_dict["top_hits_query_gff"]))

        awk_string_prefix = "awk -F'\\t' '{printf \"%s\\t%s\\t%s\\n\",$1,$4,$5}' "
        os.system(awk_string_prefix + " %s > %s" % (names_dict["top_hits_query_gff"], names_dict["top_hits_simple"]))

        stat_string = ""

        stat_string += "Reference proteins:\t%i\n" % len(reference_protein_dict)
        stat_string += "Full length top hits:\t%i\n" % full_length_top_counter
        stat_string += "Other top hits:\t%i\n" % other_top_counter
        stat_string += "Full length secondary hits:\t%i\n" % full_length_secondary_counter
        stat_string += "Other secondary hits:\t%i\n" % other_secondary_counter
        stat_string += "Proteins without hits:\t%i\n" % (len(reference_protein_dict) - full_length_top_counter - other_top_counter)

        print(stat_string)

        fd_dict["stats_general"].write(stat_string)

        for entry in fd_dict:
            fd_dict[entry].close()

    @staticmethod
    def extract_top_hits_from_target_gff(list_of_target_gff, top_hits_gff, secondary_hits_gff, id_white_list_file=None,
                                         max_hits_per_query=None):
        if id_white_list_file:
            white_ids = IdList()
            white_ids.read(id_white_list_file)
        top_hits_gff_fd = open(top_hits_gff, "w")
        secondary_hits_gff_fd = open(secondary_hits_gff, "w")
        targets_list = []
        hit_counter = 0
        gene_counter = 0
        for filename in list_of_target_gff:
            index = 0
            with open(filename, "r") as in_fd:
                #print u
                #tmp = None
                for line in in_fd:
                    tmp = line
                    if tmp == "# --- START OF GFF DUMP ---\n":
                        # read until string with target_name will appear
                        while tmp[0] == "#":
                            tmp = next(in_fd, "")

                        target_name = tmp.split("\t")[8].split(";")[1].split()[1]
                        if id_white_list_file:
                            if target_name not in white_ids:
                                continue
                        if target_name not in targets_list:
                            writing_fd = top_hits_gff_fd
                            targets_list.append(target_name)
                            gene_counter += 1
                            hit_counter = 0
                        else:
                            writing_fd = secondary_hits_gff_fd
                        # print target_name
                        hit_counter += 1
                        tmp = tmp.replace("gene_id 0", "gene_id g%i_h%i" % (gene_counter, hit_counter))
                        if hit_counter <= max_hits_per_query:
                            writing_fd.write(tmp)

                        while True:
                            tmp = next(in_fd, "")
                            # print("cccc")

                            if tmp == "# --- END OF GFF DUMP ---\n":
                                break
                            if max_hits_per_query:
                                if hit_counter > max_hits_per_query:
                                    #print "aaaaa"
                                    continue
                            writing_fd.write(tmp)
                    if tmp == "":
                        break
        top_hits_gff_fd.close()
        secondary_hits_gff_fd.close()

    @staticmethod
    def add_len_to_simple_output(top_hits_simple, len_file, out_file):
        len_dict = SynDict()
        len_dict.read(len_file)
        with open(top_hits_simple, "r") as in_fd:
            with open(out_file, "w") as out_fd:
                for line in in_fd:
                    tmp_list = line.strip().split("\t")
                    out_fd.write("%s\t%s\t%s\t%s\t%s\t%f\n" % (tmp_list[0], len_dict[tmp_list[0]], tmp_list[3],
                                                               tmp_list[1], tmp_list[2],
                                                               (float(tmp_list[2]) - float(tmp_list[1]) + 1) / float(len_dict[tmp_list[0]])))

    @staticmethod
    def extract_annotation_by_refence_id(list_of_target_gff, id_file, extracted_gff, filtered_out_gff):
        ids = IdList()
        ids.read(id_file)
        extracted_gff_fd = open(extracted_gff, "w")
        filtered_out_gff_fd = open(filtered_out_gff, "w")
        for filename in list_of_target_gff:
            with open(filename, "r") as in_fd:
                for line in in_fd:
                    tmp = line
                    if tmp == "# --- START OF GFF DUMP ---\n":
                        # read until string with target_name will appear
                        while tmp[0] == "#":
                            tmp = next(in_fd, "")

                        target_name = tmp.split("\t")[8].split(";")[1].split()[1]
                        if target_name not in ids:
                            writing_fd = filtered_out_gff_fd

                        else:
                            writing_fd = extracted_gff_fd
                        # print target_name
                        writing_fd.write(tmp)
                        while True:
                            tmp = next(in_fd, "")
                            if tmp == "# --- END OF GFF DUMP ---\n":
                                break
                            writing_fd.write(tmp)
                    if tmp == "":
                        break
        extracted_gff_fd.close()
        filtered_out_gff_fd.close()

    def prepare_annotation_file_from_transcript_and_cds(self, transcript_file, cds_file, correspondence_file,
                                                        output_prefix, format="fasta",
                                                        correspondence_key_column=0, correspondence_value_column=1,
                                                        verbose=False):
        transcript_dict = self.parse_seq_file(transcript_file, "parse", format=format)

        cds_dict = self.parse_seq_file(cds_file, "parse", format=format)

        correspondence_dict = SynDict(filename=correspondence_file, comments_prefix="#",
                                      key_index=correspondence_key_column, value_index=correspondence_value_column)

        no_corresponding_cds_transcript_list = IdList()
        cds_not_found_transcript_list = IdList()

        annotation_file = "%s.annotation" % output_prefix
        no_corresponding_cds_transcript_file = "%s.no_cds.id" % output_prefix
        cds_not_found_transcript_file = "%s.not_found_cds.id" % output_prefix

        with open(annotation_file, "w") as annotation_fd:
            for transcript_id in transcript_dict:
                if transcript_id not in correspondence_dict:
                    no_corresponding_cds_transcript_list.append(transcript_id)
                    if verbose:
                        print("No cds in correspondence file for transcript %s" % transcript_id)
                    continue
                cds_id = correspondence_dict[transcript_id]
                length = len(cds_dict[cds_id].seq)
                start = transcript_dict[transcript_id].seq.upper().find(cds_dict[cds_id].seq.upper())
                if start == -1:
                    cds_not_found_transcript_list.append(transcript_id)
                    if verbose:
                        print("CDS was not found for transcript %s" % transcript_id)
                    continue
                annotation_string = "%s\t+\t%i\t%i\n" % (transcript_id, start + 1, length)

                annotation_fd.write(annotation_string)

        no_corresponding_cds_transcript_list.write(no_corresponding_cds_transcript_file)
        cds_not_found_transcript_list.write(cds_not_found_transcript_file)

"""
    @staticmethod
    def remove_utrs_from_target_output(input_file, output_file):
        with open(input_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:
                for line in in_fd:
                    line_list = line.split("\t")
                    type = line_list[2]
                    start = int(line_list[3])
                    stop = int(line_list[4])
"""

if __name__ == "__main__":
    pass
