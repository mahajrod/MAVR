#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import SynDict, IdList


class Emapper(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "emapper.py", path=path, max_threads=max_threads)

    def parse_options(self, input, database, output_prefix, db_type=None, query_type=None, target_orthologs=None,
                      resume=None, report_orthologs=True, override_last_run=None, no_comments_in_output=None,
                      no_refine=None, no_annotation=None, no_search=None, usemem=True):

        options = " -i %s" % input
        options += " --cpu %i" % self.threads
        options += " --usemem" if usemem else ""
        options += " --temp_dir %s" % self.tmp_dir if self.tmp_dir else ""
        options += " --database %s" % database
        options += " --dbtype %s" % db_type if db_type else ""
        options += " --qtype %s" % query_type if query_type else ""
        options += " --target_orthologs %s" % target_orthologs if target_orthologs else ""

        options += " -o %s" % output_prefix if output_prefix else ""

        options += " --resume" if resume else ""
        options += " --override" if override_last_run else ""
        options += " --no_refine" if no_refine else ""
        options += " --no_annot" if no_annotation else ""
        options += " --no_search" if no_search else ""
        options += " --no_file_comments" if no_comments_in_output else ""
        options += " --report_orthologs" if report_orthologs else ""

        return options

    def assign_orthologs(self, input, database, output_prefix, db_type=None, query_type=None, target_orthologs=None,
                         resume=None, report_orthologs=True, override_last_run=None, no_comments_in_output=None,
                         no_refine=None, no_annotation=None, no_search=None, usemem=True,
                         eggnogdb_prefix=None, species_name=None, label_separator="."):

        options = self.parse_options(input, database, output_prefix,
                                     db_type=db_type,
                                     query_type=query_type,
                                     target_orthologs=target_orthologs,
                                     resume=resume,
                                     report_orthologs=report_orthologs,
                                     override_last_run=override_last_run,
                                     no_comments_in_output=no_comments_in_output,
                                     no_refine=no_refine,
                                     no_annotation=no_annotation,
                                     no_search=no_search,
                                     usemem=usemem)

        self.execute(options)

        emapper_annotation_file = "%s.emapper.annotations" % output_prefix

        self.convert_emapper_annotation_file(emapper_annotation_file,
                                             output_prefix,
                                             eggnogdb_prefix=eggnogdb_prefix,
                                             species_name=species_name,
                                             label_separator=label_separator)

    @staticmethod
    def convert_emapper_annotation_file_to_fam(emapper_annotation_file, output_fam, eggnogdb_prefix=None,
                                               species_name=None, label_separator="."):
        fam_dict = SynDict()
        with open(emapper_annotation_file, "r") as annotations_fd:
            for line in annotations_fd:
                if line[0] == "#":
                    continue
                line_list = line.split("\t")

                fam_id = line_list[10].split("|")[0]
                if not(eggnogdb_prefix is None):
                    fam_id = eggnogdb_prefix + fam_id

                gene_id = "%s%s%s" % (species_name, label_separator, line_list[0]) if species_name else line_list[0]

                if fam_id in fam_dict:
                    fam_dict[fam_id].append(gene_id)
                else:
                    fam_dict[fam_id] = [gene_id]

        fam_dict.write(filename=output_fam, splited_values=True)

    @staticmethod
    def extract_emapper_annotations_by_protein_ids(emapper_annotation_file, protein_id_file, output_annotations):
        protein_ids = IdList(filename=protein_id_file)
        with open(emapper_annotation_file, "r") as ann_fd:
            with open(output_annotations, "w") as out_fd:
                for line in ann_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    if line.split("\t")[0] in protein_ids:
                        out_fd.write(line)
    @staticmethod
    def extract_GO_terms_from_emapper_annotation_file(emapper_annotation_file, output_file):
        GO_terms_dict = SynDict(filename=emapper_annotation_file, key_index=0, value_index=5,
                                split_values=True, values_separator=",", comments_prefix="#",
                                separator="\t")
        GO_terms_dict.header = "#protein_id\tGO_terms"
        GO_terms_dict.write(output_file, header=True, splited_values=True)

        return GO_terms_dict

    @staticmethod
    def extract_predicted_gene_names_from_emapper_annotation_file(emapper_annotation_file, output_file):
        extract_predicted_gene_names_dict = SynDict(filename=emapper_annotation_file, key_index=0, value_index=4,
                                                    split_values=True, values_separator=",", comments_prefix="#",
                                                    separator="\t")
        extract_predicted_gene_names_dict.header = "#protein_id\tpredicted_gene_name"
        extract_predicted_gene_names_dict.write(output_file, header=True, splited_values=True)

        return extract_predicted_gene_names_dict

    def convert_emapper_annotation_file(self, emapper_annotation_file, output_prefix, eggnogdb_prefix=None,
                                         species_name=None, label_separator="."):
        fam_file = "%s.fam" % output_prefix
        GO_terms_file = "%s.GO" % output_prefix
        predicted_gene_names_file = "%s.predicted_gene_names" % output_prefix

        self.convert_emapper_annotation_file_to_fam(emapper_annotation_file, fam_file,
                                                    eggnogdb_prefix=eggnogdb_prefix,
                                                    species_name=species_name,
                                                    label_separator=label_separator)
        self.extract_GO_terms_from_emapper_annotation_file(emapper_annotation_file, GO_terms_file)

        self. extract_predicted_gene_names_from_emapper_annotation_file(emapper_annotation_file,
                                                                        predicted_gene_names_file)

if __name__ == "__main__":
    pass
