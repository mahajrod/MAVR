#!/usr/bin/env python
import os
import re
import time
import xmltodict

from collections import Iterable

import numpy as np

from collections import OrderedDict
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Routines.Sequence import SequenceRoutines
from CustomCollections.GeneralCollections import IdList, SynDict, TwoLvlDict, IdSet

from urllib2 import URLError


class EnsemblRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)

    @staticmethod
    def convert_biomart_protein_annotation_to_gff(biomart_annotations, output_gff, separator="\t",
                                                  extraction_mode="pfam"):

        in_fd = open(biomart_annotations, "r")
        protein_id_header = "Protein ID"

        if extraction_mode == "pfam":
            domain_id_header = "Pfam ID"
            domain_start_header = "Pfam start"
            domain_end_header = "Pfam end"

        header_list = in_fd.readline().strip().split(separator)

        for header_element in  protein_id_header, domain_id_header, domain_start_header, domain_end_header:
            if header_element not in header_list:
                raise ValueError("'%s' is absent in header of the file" % header_element)

        protein_id_index = header_list.index(protein_id_header)
        domain_id_index = header_list.index(domain_id_header)
        domain_start_index = header_list.index(domain_start_header)
        domain_end_index = header_list.index(domain_end_header)

        line_number = 1
        with open(output_gff, "w") as out_fd:
            for line in in_fd:
                line_number += 1
                tmp_list = line.strip("\n").split(separator)
                protein_id = tmp_list[protein_id_index]
                domain_id = tmp_list[domain_id_index]
                domain_start = tmp_list[domain_start_index]
                domain_end = tmp_list[domain_end_index]

                if protein_id == "":
                    print("WARNING!!! Line %i: malformed, no protein id. Skipping..." % line_number)
                    continue
                elif domain_id == "":
                    print("WARNING!!! Line %i: no domain for %s. Skipping..." % (line_number, protein_id))
                    continue
                elif domain_start == "":
                    print("WARNING!!! Line %i: no start coordinate of domain %s for protein %s. Skipping..." % (line_number, domain_id, protein_id))
                    continue
                elif domain_end == "":
                    print("WARNING!!! Line %i: no end coordinate of domain %s for protein %s. Skipping..." % (line_number, domain_id, protein_id))
                    continue

                out_fd.write("%s\tensembl\tdomain\t%s\t%s\t.\t+\t.\tID=%s; Description=%s\n" % (protein_id,
                                                                                                domain_start,
                                                                                                domain_end,
                                                                                                domain_id,
                                                                                                domain_id))
        in_fd.close()

    @staticmethod
    def get_gene_transcript_protein_from_ensembl_pep_fasta(ensembl_fasta, output_file):

        with open(ensembl_fasta, "r") as in_fd:
            with open(output_file) as out_fd:
                for line in in_fd:
                    if line[0] != ">":
                        continue

                    line_list = line.strip().split()
                    protein_id = line_list[0][1:]

                    for entry in line_list:
                        if "gene:" in entry:
                            gene_id = entry.split(":")[1]
                        elif "transcript:" in entry:
                            transcript_id = entry.split(":")[1]

                    out_fd.write("%s\t%s\t%s\n" % (gene_id, transcript_id, protein_id))

    @staticmethod
    def get_longest_pep_per_gene_from_ensembl_pep_dict(protein_dict, output_prefix=None):
        length_file = "%s.protein_length.tsv" % output_prefix
        if output_prefix:
            longest_protein_id_file = "%s.longest_pep.ids" % output_prefix

            len_fd = open(length_file, 'w')
            len_fd.write("#gene_id\tprotein_id\tprotein_length\n")

        data_dict = OrderedDict()
        for protein_id in protein_dict:
            length = protein_dict[protein_id].seq
            description_list = protein_dict[protein_id].description.split()
            print description_list
            if output_prefix:
                len_fd.write("%s\t%s\t%i\n" % (gene_id, protein_id, length))
            for entry in description_list:
                if "gene:" in entry:
                    gene_id = entry.split(":")[1]
            if gene_id not in data_dict:
                data_dict[gene_id] = protein_id
            else:
                if length > len(protein_dict[data_dict[gene_id]].seq):
                    data_dict[gene_id] = protein_id

        longest_pep_ids = IdList(data_dict.values())
        if output_prefix:
            longest_pep_ids.write(longest_protein_id_file)
            len_fd.close()
        return longest_pep_ids

    def get_longest_pep_per_gene_from_ensembl_pep_file(self, protein_file, output_prefix):
        protein_dict = self.parse_seq_file(protein_file, "parse")
        longest_pep_ids = self. get_longest_pep_per_gene_from_ensembl_pep_dict(protein_dict,
                                                                               output_prefix=output_prefix)

        longest_protein_pep_file = "%s.longest_pep.pep" % output_prefix
        SeqIO.write(self.record_by_id_generator(protein_dict, longest_pep_ids, verbose=True),
                    longest_protein_pep_file,
                    format='fasta')

        return longest_pep_ids

