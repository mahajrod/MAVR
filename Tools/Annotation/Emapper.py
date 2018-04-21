#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import SynDict, IdList


class Emapper(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "emapper.py", path=path, max_threads=max_threads)

    """
    def parse_options(self, input_fasta, kingdom, gff_file, log_file, length_cutoff=None, reject_cutoff=None,
                      evalue_cutoff=None):

        if kingdom not in self.kingdoms:
            raise ValueError("Wrong kingdom of life. Allowed: euk, bac, mito, arc")

        options = " --threads %i" % self.threads
        options += " --kingdom %s" % kingdom
        options += " --lencutoff %f" % length_cutoff if length_cutoff else ""
        options += " --reject %f" % reject_cutoff if reject_cutoff else ""
        options += " --evalue %f" % evalue_cutoff if evalue_cutoff else ""
        options += " %s" % input_fasta
        options += " > %s " % gff_file
        options += " 2> %s" % log_file

        return options
    """

    @staticmethod
    def convert_egemapper_annotation_file_to_fam(emapper_annotation_file, output_fam, eggnogdb_prefix=None,
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
    def extract_eggnogmapper_annotations_by_protein_ids(emapper_annotation_file, protein_id_file, output_annotations):
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
                                split_values=True, values_separator=True)
        GO_terms_dict.header = "#protein_id\tGO_terms"
        GO_terms_dict.write(output_file, header=True, splited_values=True)

        return GO_terms_dict

    @staticmethod
    def extract_predicted_gene_names_from_emapper_annotation_file(emapper_annotation_file, output_file):
        extract_predicted_gene_names_dict = SynDict(filename=emapper_annotation_file, key_index=0, value_index=5,
                                                    split_values=True, values_separator=True)
        extract_predicted_gene_names_dict.header = "#protein_id\tpredicted_gene_name"
        extract_predicted_gene_names_dict.write(output_file, header=True, splited_values=True)

        return extract_predicted_gene_names_dict

    def converrt_emapper_annotation_file(self, emapper_annotation_file, output_prefix, eggnogdb_prefix=None,
                                         species_name=None, label_separator="."):
        fam_file = "%s.fam" % output_prefix
        GO_terms_file = "%s.GO" % output_prefix
        predicted_gene_names_file = "%s.GO" % output_prefix

        self.convert_egemapper_annotation_file_to_fam(emapper_annotation_file, fam_file,
                                                      eggnogdb_prefix=eggnogdb_prefix,
                                                      species_name=species_name,
                                                      label_separator=label_separator)
        self.extract_GO_terms_from_emapper_annotation_file(emapper_annotation_file, GO_terms_file)

        self. extract_predicted_gene_names_from_emapper_annotation_file(emapper_annotation_file,
                                                                        predicted_gene_names_file)

if __name__ == "__main__":
    pass
