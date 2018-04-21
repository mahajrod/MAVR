#!/usr/bin/env python
import os
from copy import deepcopy
from collections import OrderedDict

from Bio import SeqIO

#from Routines import FileRoutines
#from Routines.File import FileRoutines
from Routines.Sequence import SequenceRoutines
from CustomCollections.GeneralCollections import SynDict, IdSet, IdList


class SequenceClusterRoutines(SequenceRoutines):

    def __init__(self):
        SequenceRoutines.__init__(self)

    def read_cluster_files_from_dir(self, dir_with_cluster_files):
        cluster_files_list = sorted(os.listdir(dir_with_cluster_files))
        clusters_dict = OrderedDict()
        for filename in cluster_files_list:
            filepath = "%s%s" % (self.check_path(dir_with_cluster_files), filename)
            filename_list = self.split_filename(filepath)
            clusters_dict[filename_list[1]] = SynDict()
            clusters_dict[filename_list[1]].read(filepath, header=False, separator="\t", allow_repeats_of_key=False,
                                                 split_values=True, values_separator=",", key_index=0, value_index=1,
                                                 comments_prefix="#")
        return clusters_dict

    @staticmethod
    def get_cluster_names(clusters_dict, out_file=None, white_list_ids=None):
        cluster_names = IdSet()
        for species in clusters_dict:
            species_clusters = IdSet(clusters_dict[species].keys())
            cluster_names |= species_clusters
        if out_file:
            cluster_names.write(out_file)

        return cluster_names & IdSet(white_list_ids) if white_list_ids else cluster_names

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
        cluster_names = self.get_cluster_names(clusters_dict, white_list_ids=white_list_ids)

        sequence_super_dict = OrderedDict()
        out_dir = self.check_path(output_dir)

        for species in clusters_dict:
            idx_file = "%s_tmp.idx" % species
            sequence_file = "%s%s.%s" % (self.check_path(dir_with_sequence_files), species,
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
                    #print species, cluster_id
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

        for species in clusters_dict:
            os.remove("%s_tmp.idx" % species)

    @staticmethod
    def rename_elements_in_clusters(clusters_file, syn_file, output_clusters_file,
                                    remove_clusters_with_not_renamed_elements=False,
                                    elements_with_absent_synonyms_file=None,
                                    syn_file_key_column_index=0,
                                    syn_file_value_column_index=1,
                                    syn_file_column_separator='\t'):
        syn_dict = SynDict()
        syn_dict.read(syn_file, comments_prefix="#", key_index=syn_file_key_column_index,
                      value_index=syn_file_value_column_index, separator=syn_file_column_separator)

        clusters_dict = SynDict()
        clusters_dict.read(clusters_file, split_values=True, values_separator=",", comments_prefix="#")

        output_clusters_dict = SynDict()

        absent_elements_dict = SynDict()

        for cluster in clusters_dict:
            renamed_element_list = []
            all_elements_were_renamed_flag = True
            for element in clusters_dict[cluster]:
                if element in syn_dict:
                    renamed_element_list.append(syn_dict[element])
                else:
                    if cluster not in absent_elements_dict:
                        absent_elements_dict[cluster] = [element]
                    else:
                        absent_elements_dict[cluster].append(element)
                    all_elements_were_renamed_flag = False
                    renamed_element_list.append(element)

            if (not remove_clusters_with_not_renamed_elements) or (remove_clusters_with_not_renamed_elements and all_elements_were_renamed_flag):
                output_clusters_dict[cluster] = renamed_element_list

        output_clusters_dict.write(output_clusters_file, splited_values=True)

        if elements_with_absent_synonyms_file:
            absent_elements_dict.write(elements_with_absent_synonyms_file, splited_values=True)

        return absent_elements_dict

    @staticmethod
    def check_absence_of_cluster_elements(fam_list, sequence_dict):
        absent_elements = []
        for element in fam_list:
            if element not in sequence_dict:
                absent_elements.append(element)
        return absent_elements

    def extract_sequences_from_selected_clusters(self, clusters_id_file, cluster_file, seq_file,
                                                 output_dir="./", seq_format="fasta",
                                                 out_prefix=None, create_dir_for_each_cluster=False,
                                                 skip_cluster_if_no_sequence_for_element=True,
                                                 parsing_mode="parse"):
        from Routines import SequenceRoutines
        cluster_id_list = IdList()
        cluster_dict = SynDict()
        #print(pep_file)
        self.safe_mkdir(output_dir)
        out_dir = self.check_path(output_dir)
        create_directory_for_each_cluster = True if out_prefix else create_dir_for_each_cluster
        if clusters_id_file:
            cluster_id_list.read(clusters_id_file)
        cluster_dict.read(cluster_file, split_values=True, values_separator=",")

        protein_files = self.make_list_of_path_to_files(seq_file)

        protein_dict = self.parse_seq_file(protein_files, "index_db", format=seq_format, index_file="tmp.idx") if len(protein_files) > 1 else self.parse_seq_file(protein_files[0], parsing_mode, format=seq_format, index_file="tmp.idx") #SeqIO.index_db("tmp.idx", self.make_list_of_path_to_files(seq_file), format=seq_format)

        number_of_skipped_clusters = 0
        for fam_id in cluster_id_list if clusters_id_file else cluster_dict:

            if skip_cluster_if_no_sequence_for_element:
                absent_elements = self.check_absence_of_cluster_elements(cluster_dict[fam_id], protein_dict)
                if absent_elements:
                    print "Skipping cluster %s due to absent element(%s)" % (fam_id, ",".join(absent_elements))
                    number_of_skipped_clusters += 1
                    continue

            if fam_id in cluster_dict:
                if create_directory_for_each_cluster:
                    fam_dir = "%s%s/" % (out_dir, fam_id)
                    self.safe_mkdir(fam_dir)
                    out_file = "%s%s.fasta" % (fam_dir, out_prefix if out_prefix else fam_id)
                else:
                    out_file = "%s/%s.fasta" % (out_dir, out_prefix if out_prefix else fam_id)

                SeqIO.write(SequenceRoutines.record_by_id_generator(protein_dict, cluster_dict[fam_id], verbose=True),
                            out_file, format=seq_format)

        if (len(protein_files) > 1) or (parsing_mode == "index_db"):
            os.remove("tmp.idx")
        print "%i of %i clusters were skipped due to absent elements" % (number_of_skipped_clusters, len(cluster_dict))

        return number_of_skipped_clusters

    @staticmethod
    def extract_clusters_by_element_ids(cluster_dict, element_id_list, mode="w"):
        """"
        mode: "w" - if elements from element_id_list are present in cluster extracts only that elements
              "a" - if elements from element_id_list are present in cluster extracts all elements
        """

        extracted_clusters = SynDict()
        for cluster in cluster_dict:
            extracted_elements = []
            if mode == "w":
                for element in cluster_dict[cluster]:
                    if element in element_id_list:
                        extracted_elements.append(element)
                if extracted_elements:
                    extracted_clusters[cluster] = extracted_elements
            elif mode == "a":
                for element in cluster_dict[cluster]:
                    if element in element_id_list:
                        extracted_clusters[cluster] = cluster_dict[cluster]
                        break

        return extracted_clusters

    def extract_clusters_by_element_ids_from_file(self, cluster_file, element_file, output_file, mode="w",
                                                  cluster_column=0, element_column=1, column_separator="\t",
                                                  element_separator=","):
        """"
        mode: "w" - if elements from element_id_list are present in cluster extracts only that elements
              "a" - if elements from element_id_list are present in cluster extracts all elements
        """
        cluster_dict = SynDict(filename=cluster_file, split_values=True, comments_prefix="#",
                               key_index=cluster_column, value_index=element_column,
                               separator=column_separator, values_separator=element_separator)

        element_id_list = IdList()
        element_id_list.read(element_file, comments_prefix="#")
        extracted_clusters = self.extract_clusters_by_element_ids(cluster_dict, element_id_list, mode=mode)
        extracted_clusters.write(output_file, splited_values=True)

    @staticmethod
    def label_cluster_elements(cluster_dict, label, separator="@", label_position="first"):
        labeled_cluster_dict = SynDict()
        if label_position == "first":
            label_function = lambda s: "%s%s%s" % (label, separator, s)
        elif label_position == "last":
            label_function = lambda s: "%s%s%s" % (s, separator, label)

        for cluster in cluster_dict:
            labeled_cluster_dict[cluster] = []
            for element in cluster_dict[cluster]:
                labeled_cluster_dict[cluster].append(label_function(element))

        return labeled_cluster_dict

    def label_cluster_elements_from_file(self, input_file, label, output_file, separator="@", label_position="first"):
        input_dict = SynDict()
        input_dict.read(input_file, split_values=True, comments_prefix="#")

        output_dict = self.label_cluster_elements(input_dict, label, separator=separator,
                                                  label_position=label_position)
        output_dict.write(output_file, splited_values=True)

        return output_dict

    @staticmethod
    def extract_labels_from_cluster_elements(cluster_dict, separator="@", label_position="first"):
        labels_set = set()
        for cluster in cluster_dict:
            for element in cluster_dict[cluster]:
                label = element.split(separator)
                label = label[0] if label_position == "first" else label[-1]
                labels_set.add(label)
        return labels_set

    @staticmethod
    def replace_label(cluster_dict, syn_dict=None, old_separator="@", old_label_position="first", new_separator="@", new_label_position="first"):
        new_cluster_dict = SynDict()
        for cluster in cluster_dict:
            new_cluster_dict[cluster] = []
            for element in cluster_dict[cluster]:
                tmp = element.split(old_separator)
                if old_label_position == "first":
                    label = tmp[0]
                    element_id = old_separator.join(tmp[1:])
                else:
                    label = tmp[-1]
                    element_id = old_separator.join(tmp[:-1])

                if new_label_position == 'first':
                    new_cluster_dict[cluster].append("%s%s%s" % (syn_dict[label] if syn_dict else label, new_separator, element_id))
                else:
                    new_cluster_dict[cluster].append("%s%s%s" % (element_id, new_separator, syn_dict[label] if syn_dict else label))

        return new_cluster_dict

    def replace_label_from_file(self, input_file, output_file, syn_file_or_dict, old_separator="@",
                                old_label_position="first",
                                new_separator="@", new_label_position="first"):

        syn_dict = SynDict(filename=syn_file_or_dict, split_values=False) if isinstance(syn_file_or_dict, str) else SynDict(syn_file_or_dict)
        cluster_dict = SynDict(filename=input_file, split_values=True)
        new_cluster_dict = self.replace_label(cluster_dict, syn_dict=syn_dict, old_separator=old_separator,
                                              old_label_position=old_label_position, new_separator=new_separator,
                                              new_label_position=new_label_position)

        new_cluster_dict.write(output_file, splited_values=True)

        return new_cluster_dict

    @staticmethod
    def extract_single_copy_clusters(dict_of_cluster_dicts, label_elements=False, separator="@",
                                     label_position="first"):

        if label_position == "first":
            label_function = lambda s, label: "%s%s%s" % (label, separator, s)
        elif label_position == "last":
            label_function = lambda s, label: "%s%s%s" % (s, separator, label)

        sc_clusters_dict = SynDict()

        clusters_set = set()
        for group in dict_of_cluster_dicts:
            clusters_set = clusters_set | set(dict_of_cluster_dicts[group].keys())

        for cluster in clusters_set:
            for group in dict_of_cluster_dicts:
                if cluster not in dict_of_cluster_dicts[group]:
                    break
                if len(dict_of_cluster_dicts[group][cluster]) > 1:
                    break
            else:
                sc_clusters_dict[cluster] = []
                for group in dict_of_cluster_dicts:
                    if label_elements:
                        sc_clusters_dict[cluster].append(label_function(dict_of_cluster_dicts[group][cluster][0],
                                                                        group))
                    else:
                        sc_clusters_dict[cluster].append(dict_of_cluster_dicts[group][cluster][0])

        return sc_clusters_dict

    def extract_single_copy_clusters_from_files(self, list_of_cluster_files, output_file, label_elements=False, separator="@",
                                                label_position="first", function_to_convert_filename_to_label=None):
        dict_of_cluster_dicts = OrderedDict()
        for filename in list_of_cluster_files:
            if function_to_convert_filename_to_label:
                label = function_to_convert_filename_to_label(filename)
            else:
                label = self.split_filename(filename)[1]  # use basename as label

            dict_of_cluster_dicts[label] = SynDict()
            dict_of_cluster_dicts[label].read(filename, split_values=True, comments_prefix="#")

        sc_clusters_dict = self. extract_single_copy_clusters(dict_of_cluster_dicts, label_elements=label_elements,
                                                              separator=separator, label_position=label_position)

        sc_clusters_dict.write(output_file, splited_values=True)

        return sc_clusters_dict

    @staticmethod
    def remove_elements_by_ids(input_dict, black_list, mode="full"):
        filtered_dict = SynDict()

        for cluster_id in input_dict:
            elements = []
            for element_id in input_dict[cluster_id]:
                if mode == "full":
                    if element_id not in black_list:
                        elements.append(element_id)
                elif mode == "partial":
                    for black_id_part in black_list:
                        #print black_list
                        if black_id_part in element_id:
                            #print black_id_part, element_id
                            break
                    else:
                        #print "aaaa"
                        elements.append(element_id)
                else:
                    raise ValueError("Unknown mode for id comparison")
            if elements:
                filtered_dict[cluster_id] = elements

        return filtered_dict

    def remove_elements_by_ids_from_files(self, input_file, output_file, black_list_file, mode="full"):

        cluster_dict = SynDict(filename=input_file, split_values=True)
        black_list = IdList(filename=black_list_file)

        filtered_dict = self.remove_elements_by_ids(cluster_dict, black_list, mode=mode)

        filtered_dict.write(output_file, splited_values=True)

    @staticmethod
    def split_element_id(element_id, separator="@", label_position="first"):
        element_list = element_id.split(separator)
        if len(element_list) > 2:
            raise ValueError("Multiple separators in element id(%s)" % element_id)
        elif len(element_list) < 2:
            raise ValueError("No separator was found in element id(%s)" % element_id)
        if label_position == "first":
            return element_list
        elif label_position == "last":
            return [element_list[1], element_list[0]]

    def extract_clusters_and_elements_by_labels(self, cluster_dict, label_list, separator="@", label_position="first"):
        filtered_dict = SynDict()

        for cluster_id in cluster_dict:
            for element_id in cluster_dict[cluster_id]:
                label, element = self.split_element_id(element_id, separator=separator,
                                                       label_position=label_position)
                if label in label_list:
                    if cluster_id in filtered_dict:
                        filtered_dict[cluster_id].append(element_id)
                    else:
                        filtered_dict[cluster_id] = [element_id]

        return filtered_dict

    def extract_clusters_and_elements_by_labels_from_files(self, cluster_file, label_file, output_file,
                                                           separator="@", label_position="first"):
        cluster_dict = SynDict(filename=cluster_file, split_values=True)
        label_list = IdList(filename=label_file) if isinstance(label_file, str) else label_file

        output_dict = self.extract_clusters_and_elements_by_labels(cluster_dict, label_list, separator=separator,
                                                                   label_position=label_position)

        output_dict.write(output_file)




