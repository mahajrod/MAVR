#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from Tools.Abstract import Tool

from Tools.LinuxTools import CGAS
from Parsers.TRF import CollectionTRF
from Routines.File import split_filename, save_mkdir
from CustomCollections.GeneralCollections import SynDict


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
                           store_intermediate_files=False,
                           splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                           number_of_results_to_report=None,
                           other_options=None,
                           num_of_files=None,
                           converted_output_dir="converted_output"):
        splited_filename = split_filename(query_file)
        self.split_fasta(query_file, splited_fasta_dir, num_of_recs_per_file=num_of_recs_per_file,
                         num_of_files=num_of_files,
                         output_prefix=splited_filename[1])

        common_options = self.parse_common_options(model, show_alignment=show_alignment,
                                                   show_sugar=show_sugar, show_cigar=show_cigar,
                                                   show_vulgar=show_vulgar, show_query_gff=show_query_gff,
                                                   show_target_gff=show_target_gff,
                                                   number_of_results_to_report=number_of_results_to_report,
                                                   other_options=other_options)

        options_list = []
        splited_files = os.listdir(splited_fasta_dir)

        save_mkdir(splited_result_dir)
        #save_mkdir(converted_output_dir)

        for filename in splited_files:
            filename_list = split_filename(filename)
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
                            fd_dict["gff"].write(tmp)
                            if index == 0:
                                fd_dict["query_gff"].write(tmp)
                            else:
                                fd_dict["target_gff"].write(tmp)
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

        awk_string_prefix = "awk -F'\t' '{printf \"%s\t%s\t%s\n\",$1,$4,$5}' "
        os.system(awk_string_prefix + " %s > %s" % (top_hits_query_gff, top_hits_simple))

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

if __name__ == "__main__":
    pass
