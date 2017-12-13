#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from collections import OrderedDict
from Tools.Abstract import Tool

from CustomCollections.GeneralCollections import SynDict, IdList


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
                             number_of_results_to_report=None, other_options=None):

        options = " --model %s" % model
        options += " --showalignment" if show_alignment else ""
        options += " --showsugar" if show_sugar else ""
        options += " --showcigar" if show_cigar else ""
        options += " --showvulgar" if show_vulgar else ""
        options += " --showquerygff" if show_query_gff else ""
        options += " --showtargetgff" if show_target_gff else ""
        options += " -Q %s" % query_type if query_type else ""
        options += " -T %s" % target_type if target_type else ""
        options += " -n %i" % number_of_results_to_report if number_of_results_to_report else ""
        options += " %s" % other_options if other_options else ""

        return options

    def parallel_alignment(self, query_file, target_file, model, num_of_recs_per_file=None,
                           show_alignment=None, show_sugar=True, show_cigar=None,
                           show_vulgar=None, show_query_gff=None, show_target_gff=None,
                           store_intermediate_files=True,
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
                                                   other_options=other_options)

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

        """
        for filename in splited_files:

            trf_output_file = "%s/%s.%i.%i.%i.%i.%i.%i.%i.dat" % (splited_result_dir, filename,
                                                                  matching_weight, mismatching_penalty,
                                                                  indel_penalty, match_probability,
                                                                  indel_probability,
                                                                  min_alignment_score, max_period)

            self.convert_trf_report(trf_output_file, "%s/%s" % (converted_output_dir, filename))

        for suffix in (".rep", ".gff", ".simple.gff", ".short.tab", ".wide.tab"):
            file_str = ""
            merged_file = "%s%s" % (output_prefix, suffix)
            for filename in splited_files:
                file_str += " %s/%s%s" % (converted_output_dir, filename, suffix)
            CGAS.cat(file_str, merged_file)
        """
        if not store_intermediate_files:
            shutil.rmtree(splited_fasta_dir)
            #shutil.rmtree(splited_result_dir)
            #shutil.rmtree(converted_output_dir)

    @staticmethod
    def split_output(exonerate_output_files, output_prefix):
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
                      }
        top_hits_gff = "%s.top_hits.gff" % output_prefix
        top_hits_vulgar = "%s.top_hits.vulgar" % output_prefix
        top_hits_sugar = "%s.top_hits.sugar" % output_prefix
        top_hits_target_gff = "%s.top_hits.target.gff" % output_prefix
        top_hits_query_gff = "%s.top_hits.query.gff" % output_prefix
        top_hits_simple = "%s.top_hits.query.simple" % output_prefix

        fd_dict = {}
        for output_type in names_dict:
            fd_dict[output_type] = open(names_dict[output_type], "w")

        for filename in exonerate_output_files:
            index = 0
            with open(filename, "r") as in_fd:
                #print u
                #tmp = None
                for line in in_fd:
                    tmp = line
                    if tmp[:13] == "C4 Alignment:":
                        tmp = in_fd.next()
                        fd_dict["alignment"].write(tmp)
                        while True:
                            tmp = next(in_fd, "")
                            if (tmp[0] == "c") or (tmp[0] == "v") or (tmp[0] == "s") or (tmp[0] == "#"):
                                break
                            fd_dict["alignment"].write(tmp)
                            if tmp == "":
                                break
                    if tmp == "# --- START OF GFF DUMP ---\n":
                        fd_dict["gff"].write(tmp)
                        if index == 0:

                            fd_dict["query_gff"].write(tmp)
                        else:
                            fd_dict["target_gff"].write(tmp)
                        while True:
                            tmp = next(in_fd, "")
                            if "\tgene\t" in tmp:
                                current_gene_id = tmp.strip().split("\t")[-1].split(";").split()[1]
                                print current_gene_id
                            elif tmp[0] != "#":
                                tmp = tmp.strip() + "; gene_id %s\n" % current_gene_id
                            fd_dict["gff"].write(tmp)
                            if index == 0:
                                fd_dict["query_gff"].write(tmp)
                            else:
                                fd_dict["target_gff"].write(tmp)
                            if tmp[0] != "#":
                                if "\tsplice" in tmp:
                                    fd_dict["splice"].write(tmp)
                                elif "\texon\t" in tmp:
                                    fd_dict["exon"].write(tmp.strip() + "; transcript_id %s\n" % current_gene_id)
                                elif "\tintron\t" in tmp:
                                    fd_dict["intron"].write(tmp)
                                elif "\tcds\t" in tmp:
                                    fd_dict["cds"].write(tmp.strip() + "; transcript_id %s\n" % current_gene_id)
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
                    elif tmp[:6] == "sugar:":
                        fd_dict["sugar"].write(tmp[6:])
                    elif tmp[:6] == "cigar:":
                        fd_dict["cigar"].write(tmp[6:])
        for output_type in fd_dict:
            fd_dict[output_type].close()

        # extract top hits from vulgar
        awk_string_prefix = "awk '{if (substr($0,0,1) != \"#\") {if ($1 != SEQ_ID) {print $0}; SEQ_ID=$1}}' "
        os.system(awk_string_prefix + " %s > %s" % (names_dict["vulgar"], top_hits_vulgar))
        os.system(awk_string_prefix + " %s > %s" % (names_dict["sugar"], top_hits_sugar))
        os.system(awk_string_prefix + " %s > %s" % (names_dict["query_gff"], top_hits_query_gff))

        awk_string_prefix = "awk -F'\\t' '{printf \"%s\\t%s\\t%s\\n\",$1,$4,$5}' "
        os.system(awk_string_prefix + " %s > %s" % (top_hits_query_gff, top_hits_simple))

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
