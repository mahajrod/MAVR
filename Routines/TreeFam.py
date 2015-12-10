#!/usr/bin/env python
import os

from Bio import SeqIO

from Routines.File import save_mkdir, check_path
from Routines.Sequence import record_by_id_generator
from CustomCollections.GeneralCollections import IdList, SynDict


class TreeFamRoutines:
    def __init__(self):
        pass

    @staticmethod
    def extract_proteins_from_selected_families(families_id_file, fam_file, pep_file,
                                                output_dir="./", pep_format="fasta",
                                                out_prefix=None, create_dir_for_each_family=False):
        fam_id_list = IdList()
        fam_dict = SynDict()

        save_mkdir(output_dir)
        out_dir = check_path(output_dir)
        create_directory_for_each_family = True if out_prefix else create_dir_for_each_family

        fam_id_list.read(families_id_file)
        fam_dict.read(fam_file, split_values=True, values_separator=",")
        protein_dict = SeqIO.index_db("tmp.idx", pep_file, format=pep_format)

        for fam_id in fam_id_list:
            if fam_id in fam_dict:
                if create_directory_for_each_family:
                    fam_dir = "%s%s/" % (out_dir, fam_id)
                    save_mkdir(fam_dir)
                    out_file = "%s%s.pep" % (fam_dir, out_prefix if out_prefix else fam_id)
                else:
                    out_file = "%s.pep" % (out_prefix if out_prefix else fam_id)

                SeqIO.write(record_by_id_generator(protein_dict, fam_dict[fam_id]), out_file, format=pep_format)

        os.remove("tmp.idx")

    @staticmethod
    def add_length_to_fam_file(fam_file, len_file, out_file):
        fam_dict = SynDict()
        fam_dict.read(fam_file, split_values=True, comments_prefix="#")
        len_dict = SynDict()
        len_dict.read(len_file, comments_prefix="#")

        with open(out_file, "w") as out_fd:
            for family in fam_dict:
                len_list = []
                for member in fam_dict[family]:
                    len_list.append(None if member not in len_dict else len_dict[member])

                out_fd.write("%s\t%s\t%s\n" % (family, ",".join(fam_dict[family]), ",".join(len_list)))

