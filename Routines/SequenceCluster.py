#!/usr/bin/env python
import os
from copy import deepcopy
from collections import OrderedDict

from Bio import SeqIO

from Routines import FileRoutines
from Routines.Sequence import SequenceRoutines
from CustomCollections.GeneralCollections import SynDict, IdSet


class SequenceClusterRoutines:
    def __init__(self):

        pass

    """
    @staticmethod
    def split_proteins_per_species(dir_with_proteins, output_dir, input_format="fasta", output_format="fasta"):
        input_files = FileRoutines.make_list_of_path_to_files([dir_with_proteins] if isinstance(dir_with_proteins, str) else dir_with_proteins)

        out_dir = FileRoutines.check_path(output_dir)
        FileRoutines.save_mkdir(out_dir)

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
    """

    @staticmethod
    def read_cluster_files_from_dir(dir_with_cluster_files):
        cluster_files_list = sorted(os.listdir(dir_with_cluster_files))
        clusters_dict = OrderedDict()
        for filename in cluster_files_list:
            filepath = "%s%s" % (FileRoutines.check_path(dir_with_cluster_files), filename)
            filename_list = FileRoutines.split_filename(filepath)
            clusters_dict[filename_list[1]] = SynDict()
            clusters_dict[filename_list[1]].read(filepath, header=False, separator="\t", allow_repeats_of_key=False,
                                                 split_values=True, values_separator=",", key_index=0, value_index=1,
                                                 comments_prefix="#")
        return clusters_dict

    @staticmethod
    def get_cluster_names(clusters_dict, out_file=None):
        cluster_names = IdSet()
        for species in clusters_dict:
            species_clusters = IdSet(clusters_dict[species].keys())
            cluster_names |= species_clusters
        if out_file:
            cluster_names.write(out_file)
        return cluster_names

    @staticmethod
    def get_sequence_names(clusters_dict, write_ids=False, out_prefix=None, white_list_ids=None):
        sequence_names_dict = SynDict()
        for species in clusters_dict:
            sequence_names_dict[species] = IdSet()
        for species in clusters_dict:
            for cluster_id in clusters_dict[species]:
                if white_list_ids:
                    if cluster_id not in white_list_ids:
                        continue
                sequence_names_dict[species] = sequence_names_dict[species] | IdSet(clusters_dict[species][cluster_id])
        if write_ids:
            for species in clusters_dict:
                out_file = "%s_%s.ids" % (out_prefix, species) if out_prefix else "%s.ids" % species
                sequence_names_dict[species].write(out_file)
        return sequence_names_dict

    @staticmethod
    def merge_clusters(clusters_dict, label_species="False", separator_for_labeling="_",
                       species_label_first=True):

        if species_label_first:
            label_sequence = lambda label, name: "%s%s%s" % (label, separator_for_labeling, name)
        else:
            label_sequence = lambda label, name: "%s%s%s" % (name, separator_for_labeling, label)
        if label_species:
            expression = label_sequence
        else:
            expression = lambda label, name: name

        merged_clusters = SynDict()
        for species in clusters_dict:
            for cluster in clusters_dict[species]:
                if cluster not in merged_clusters:
                    merged_clusters[cluster] = []
                for sequence_name in clusters_dict[species][cluster]:
                    merged_clusters[cluster].append(expression(species, sequence_name))

        return merged_clusters

    def merge_clusters_from_files(self, dir_with_cluster_files, output_file, label_species="False",
                                  separator_for_labeling="@", species_label_first=True):

        clusters_dict = self.read_cluster_files_from_dir(dir_with_cluster_files)
        merged_clusters = self.merge_clusters(clusters_dict, label_species=label_species,
                                              separator_for_labeling=separator_for_labeling,
                                              species_label_first=species_label_first)
        merged_clusters.write(output_file, splited_values=True)
        return merged_clusters

    def extract_monocluster_ids(self, clusters_dict, white_list_ids=None, out_file=None):
        """
        Extracts clusters with only one sequence in all species.
        """
        monocluster_ids = IdSet()
        cluster_names = self.get_cluster_names(clusters_dict)

        for cluster_name in cluster_names:
            for species in clusters_dict:
                if white_list_ids:
                    if cluster_name not in white_list_ids:
                        break
                if cluster_name not in clusters_dict[species]:
                    break
                if len(clusters_dict[species][cluster_name]) > 1:
                    break
            else:
                monocluster_ids.add(cluster_name)

        if out_file:
            monocluster_ids.write(out_file)

        return monocluster_ids

    def extract_monocluster_ids_from_file(self, dir_with_cluster_files, out_file, file_with_white_list_ids=None):
        # filenames are counted as species names
        white_list_ids = None
        if file_with_white_list_ids:
            white_list_ids = IdSet()
            white_list_ids.read(file_with_white_list_ids)
        clusters_dict = self.read_cluster_files_from_dir(dir_with_cluster_files)
        monoclusters = self.extract_monocluster_ids(clusters_dict, out_file=out_file, white_list_ids=white_list_ids)
        return monoclusters

    def extract_sequences_by_clusters(self, dir_with_cluster_files, dir_with_sequence_files, output_dir,
                                      file_with_white_list_cluster_ids=None, mode="families",
                                      sequence_file_extension="fasta", sequence_file_format="fasta",
                                      label_species=False, separator_for_labeling="@", species_label_first=True):
        """
        basenames of cluster and sequence files must be same

        mode:
            clusters - extract sequences from clusters in separate files,
            species - extract sequences from species to separate files
        """
        white_list_ids = None
        if file_with_white_list_cluster_ids:
            white_list_ids = IdSet()
            white_list_ids.read(file_with_white_list_cluster_ids)

        clusters_dict = self.read_cluster_files_from_dir(dir_with_cluster_files)
        cluster_names = self.get_cluster_names(clusters_dict)

        sequence_super_dict = OrderedDict()
        out_dir = FileRoutines.check_path(output_dir)

        for species in clusters_dict:
            idx_file = "%s_tmp.idx" % species
            sequence_file = "%s%s.%s" % (FileRoutines.check_path(dir_with_sequence_files), species,
                                         sequence_file_extension)
            sequence_super_dict[species] = SeqIO.index_db(idx_file, sequence_file, format=sequence_file_format)

        if mode == "species":
            seqeuence_names = self.get_sequence_names(clusters_dict, write_ids=False, out_prefix=None,
                                                      white_list_ids=white_list_ids)
            for species in seqeuence_names:
                out_file = "%s%s.%s" % (out_dir, species, sequence_file_extension)
                SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_super_dict[species],
                                                                    seqeuence_names[species]),
                            out_file, format=sequence_file_format)
        elif mode == "families":

            def per_family_record_generator(seq_super_dict, clust_dict, cluster_id):
                if species_label_first:
                    label_sequence = lambda label, name: "%s%s%s" % (label, separator_for_labeling, name)
                else:
                    label_sequence = lambda label, name: "%s%s%s" % (name, separator_for_labeling, label)

                for species in seq_super_dict:
                    for record_id in clust_dict[species][cluster_id]:
                        if label_species:
                            record = deepcopy(seq_super_dict[species][record_id])
                            record.id = label_sequence(species, record_id)
                            yield record
                        else:
                            yield seq_super_dict[species][record_id]

            for cluster_name in cluster_names:
                out_file = "%s%s.%s" % (out_dir, cluster_name, sequence_file_extension)
                SeqIO.write(per_family_record_generator(sequence_super_dict, clusters_dict, cluster_name),
                            out_file, format=sequence_file_format)





