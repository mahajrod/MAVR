#!/usr/bin/env python
import os

from copy import deepcopy

from Bio import SeqIO

from Routines import NCBIRoutines
#from Routines.File import make_list_of_path_to_files
from CustomCollections.GeneralCollections import IdList, SynDict
from Routines.SequenceCluster import SequenceClusterRoutines


class EggNOGRoutines(SequenceClusterRoutines):
    def __init__(self):
        SequenceClusterRoutines.__init__(self)

    @staticmethod
    def edit_profile_names_in_fam_file(input_fam_file, output_fam_file, comments_prefix="#"):
        if comments_prefix:
            comments_prefix_len = len(comments_prefix)

        with open(input_fam_file, "r") as in_fd:
            with open(output_fam_file, "w") as out_fd:
                for line in in_fd:
                    if comments_prefix:
                        if line[:comments_prefix_len] == comments_prefix:
                            out_fd.write(line)
                    line_list = line.split("\t")
                    fam_name = line_list[0].split(".")[1]
                    out_fd.write("%s\t%s" % (fam_name, line_list[1]))

    @staticmethod
    def convert_members_tsv_to_fam(input_file, output_file):
        cmd = "awk -F'\t' '{printf \"%%s\\t%%s\\n\",$2,$6 }' %s > %s" % (input_file, output_file)
        os.system(cmd)

    def extract_proteins_from_alignments(self, dir_with_alignments, output_dir):
        out_dir = self.check_path(output_dir)

        #print type(FileRoutines)

        input_files = self.make_list_of_path_to_files([dir_with_alignments] if isinstance(dir_with_alignments, str) else dir_with_alignments)

        self.safe_mkdir(out_dir)
        from Routines import MultipleAlignmentRoutines
        for filename in input_files:
            filename_list = self.split_filename(filename)
            output_file = "%s%s%s" % (out_dir, filename_list[1], filename_list[2])
            MultipleAlignmentRoutines.extract_sequences_from_alignment(filename, output_file)

    def split_proteins_per_species(self, dir_with_proteins, output_dir, input_format="fasta", output_format="fasta"):
        #print type(FileRoutines)
        input_files = self.make_list_of_path_to_files([dir_with_proteins] if isinstance(dir_with_proteins, str) else dir_with_proteins)

        out_dir = self.check_path(output_dir)
        self.safe_mkdir(out_dir)

        protein_dict = SeqIO.index_db("temp.idx", input_files, format=input_format)

        syn_dict = SynDict()

        for protein_id in protein_dict:
            taxa_id = protein_id.split(".")[0]
            # pep_id = ".".join(tmp_list[1:])
            if taxa_id not in syn_dict:
                syn_dict[taxa_id] = []
            syn_dict[taxa_id].append(protein_id)

        def renamed_records_generator(record_dict, taxa_id):
            for record_id in syn_dict[taxa_id]:
                record = deepcopy(record_dict[record_id])
                #print(record)
                record.id = ".".join(record_id.split(".")[1:])
                yield record

        for taxa_id in syn_dict:
            out_file = "%s%s.pep" % (out_dir, taxa_id)
            SeqIO.write(renamed_records_generator(protein_dict, taxa_id), out_file, format=output_format)

    def get_species_from_eggnog_tsv(self, eggnog_tsv, output_prefix, email=None):

        cluster_dict = SynDict(filename=eggnog_tsv, key_index=1, value_index=5, split_values=True)

        species_ids = self.extract_labels_from_cluster_elements(cluster_dict, separator=".", label_position="first")
        print("Input species ids: %s" % " ".join(species_ids))

        if not email:
            species = species_ids
        else:
            species = NCBIRoutines.get_taxonomy(species_ids, "%s.species.taxonomy" % output_prefix,
                                                email, input_type="id")

        species.write("%s.species" % output_prefix, splited_values=True)

        for species_id in species:
            for i in range(0, len(species[species_id])):
                species[species_id][i] = species[species_id][i].lower().replace(" ", "_")

        species.write("%s.replaced_spaces.species" % output_prefix, splited_values=True)

    def extract_eggnog_fam_by_protein_syn_dict(self, eggnog_fam_dict, protein_syn_dict, output_prefix=None, species_id=None):

        extracted_families = SynDict()
        common_protein_names_to_families_dict = SynDict()
        common_names_to_eggnog_proteins_syn_dict = SynDict()

        not_found_proteins_common_names = IdList()

        transposed_eggnog_fam_dict = eggnog_fam_dict.exchange_key_and_value()

        for common_protein_name in protein_syn_dict:
            not_found = True
            for protein_id in protein_syn_dict[common_protein_name]:
                extended_protein_id = protein_id if species_id is None else species_id + "." + protein_id
                if extended_protein_id in transposed_eggnog_fam_dict:
                    not_found = False
                    if common_protein_name not in common_protein_names_to_families_dict:
                        common_protein_names_to_families_dict[common_protein_name] = [transposed_eggnog_fam_dict[extended_protein_id][0]]
                        common_names_to_eggnog_proteins_syn_dict[common_protein_name] = [extended_protein_id]
                    else:
                        common_protein_names_to_families_dict[common_protein_name].append(transposed_eggnog_fam_dict[extended_protein_id][0])
                        common_names_to_eggnog_proteins_syn_dict[common_protein_name].append(extended_protein_id)
                    if transposed_eggnog_fam_dict[extended_protein_id][0] not in extracted_families:
                        extracted_families[transposed_eggnog_fam_dict[extended_protein_id][0]] = eggnog_fam_dict[transposed_eggnog_fam_dict[extended_protein_id][0]]

            if not_found:
                not_found_proteins_common_names.append(common_protein_name)

        if output_prefix:
            extracted_families.write(filename="%s.extracted_families.fam" % output_prefix, splited_values=True)
            common_protein_names_to_families_dict.write(filename="%s.common_protein_names_to_families.correspondence" % output_prefix, splited_values=True)
            common_names_to_eggnog_proteins_syn_dict.write(filename="%s.common_protein_names_to_eggnog_proteins.correspondence" % output_prefix, splited_values=True)
            not_found_proteins_common_names.write(filename="%s.not_found.common_names" % output_prefix)

            #print common_names_to_eggnog_proteins_syn_dict
            #print common_protein_names_to_families_dict
        return extracted_families, common_protein_names_to_families_dict, \
               common_names_to_eggnog_proteins_syn_dict, not_found_proteins_common_names

