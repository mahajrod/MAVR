#!/usr/bin/env python
import os

from copy import deepcopy

from Bio import SeqIO

from Routines import FileRoutines, NCBIRoutines
from Routines.File import make_list_of_path_to_files
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

    @staticmethod
    def extract_proteins_from_alignments(dir_with_alignments, output_dir):
        out_dir = FileRoutines.check_path(output_dir)

        print type(FileRoutines)

        input_files = make_list_of_path_to_files([dir_with_alignments] if isinstance(dir_with_alignments, str) else dir_with_alignments)

        FileRoutines.safe_mkdir(out_dir)
        from Routines import MultipleAlignmentRoutines
        for filename in input_files:
            filename_list = FileRoutines.split_filename(filename)
            output_file = "%s%s%s" % (out_dir, filename_list[1], filename_list[2])
            MultipleAlignmentRoutines.extract_sequences_from_alignment(filename, output_file)

    @staticmethod
    def split_proteins_per_species(dir_with_proteins, output_dir, input_format="fasta", output_format="fasta"):
        input_files = FileRoutines.make_list_of_path_to_files([dir_with_proteins] if isinstance(dir_with_proteins, str) else dir_with_proteins)

        out_dir = FileRoutines.check_path(output_dir)
        FileRoutines.safe_mkdir(out_dir)

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
        print species_ids
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








