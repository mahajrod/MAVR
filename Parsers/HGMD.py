#!/usr/bin/env python
__author__ = 'mahajrod'

import numpy as np
import pandas as pd


class HGMDvariants(pd.DataFrame):
    def __init__(self, input, separator="\t", header=None):

        pd.DataFrame.__init__(self)

        if isinstance(input, str) or isinstance(input, file):
            pd.DataFrame.__init__(self, pd.read_csv(input, sep=separator, header=0 if header is True else header))
        elif input:
            pd.DataFrame.__init__(self, input)
        self.disease_column_name = "disease"
        self.gene_column_name = "gene"
        self.hgvs_column_name = "hgvsall"
        self.codon_change_column = "base"
        self.aa_change_column = "amino"
        self.changed_codon_position = "codon"

        self.hgvs_snp_template1 = '(\d+)([A-Za-z])to([A-Za-z])'
        self.hgvs_snp_template2 = '([A-Za-z])(\d+)([A-Za-z*])'

        self.hgvs_indel_template1 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(del|dup|ins)([A-Za-z]+)'
        self.hgvs_indel_template2 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(del|dup|ins)(\d+)'
        self.hgvs_indel_template3 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(del|dup|ins)'

        self.hgvs_hgmd_indel_template1 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(del)([A-Za-z]+)(ins)([A-Za-z]+)?'
        self.hgvs_hgmd_indel_template2 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(del)(\d+)(ins)(\d+)?'
        self.hgvs_hgmd_indel_template3 = '(-?\d+[\-\+]?\d*)_?(-?\d+[\-\+]?\d*)?(delins)'

    def get_snp_hgvs_dataframes(self):
        tmp = self[self.hgvs_column_name].str.split("|", expand=True)
        return tmp[0].str.extract(self.hgvs_snp_template1), \
               tmp[1].str.extract(self.hgvs_snp_template2)

    def get_codon_change_dataframe(self):
        return self[self.codon_change_column].str.split("-", expand=True)

    def get_aa_change_dataframe(self):
        return self[self.aa_change_column].str.split("-", expand=True)

    def get_changed_codon_position_dataframe(self):
        return self[self.changed_codon_position]

    def get_indel_hgvs_dataframes(self):
        tmp = self[self.hgvs_column_name].str.split("|", expand=True)
        return tmp[0].str.extract(self.hgvs_indel_template1), \
               tmp[1].str.extract(self.hgvs_indel_template2), \
               tmp[2].str.extract(self.hgvs_indel_template3)

    def get_hgmd_indel_hgvs_dataframes(self):
        tmp = self[self.hgvs_column_name].str.split("|", expand=True)
        return tmp[0].str.extract(self.hgvs_hgmd_indel_template1), \
               tmp[1].str.extract(self.hgvs_hgmd_indel_template2), \
               tmp[2].str.extract(self.hgvs_hgmd_indel_template3)

    def extract_snp_coordinates(self):
        hgvs_nuc, hgvs_aa = self.get_snp_hgvs_dataframes()

        results_df = pd.concat([self[self.disease_column_name],               # 1 column
                                self[self.gene_column_name],                  # 1 column
                                hgvs_nuc[0],                                  # 1 column
                                hgvs_nuc[1],                                  # 1 column
                                hgvs_nuc[2],                                  # 1 column
                                self.get_codon_change_dataframe(),            # 2 columns
                                self.get_changed_codon_position_dataframe(),  # 1 column
                                hgvs_aa[0],                                   # 1 column
                                hgvs_aa[2],                                   # 1 column
                                self.get_aa_change_dataframe()                # 2 columns
                                ],
                               axis=1
                               )

        results_df.columns = ["disease",
                              "gene",
                              "nuc_pos",
                              "ref_allel",
                              "alt_allel",
                              "ref_codon",
                              "alt_codon",
                              "codon_pos",
                              "ref_aa,1l",
                              "alt_aa,1l",
                              "ref_aa,3l",
                              "alt_aa,3l"]

        return results_df

    def extract_indel_coordinates(self):
        hgvs_indel_1, hgvs_indel_2, hgvs_indel_3 = self.get_indel_hgvs_dataframes()

        results_df = pd.concat([self[self.disease_column_name],                 # 1 column
                                self[self.gene_column_name],                    # 1 column
                                hgvs_indel_1[2],                                # 1 column
                                hgvs_indel_1[0],                                # 1 column
                                hgvs_indel_1[1],                                # 1 column
                                hgvs_indel_1[3],                                # 1 column
                                hgvs_indel_2[3],                                # 1 column
                                self.get_changed_codon_position_dataframe(),    # 1 column
                                ],
                               axis=1
                               )

        results_df.columns = ["disease",
                              "gene",
                              "indel_type",
                              "indel_start",
                              "indel_end",
                              "indel_sequence",
                              "indel_length",
                              "codon_pos"]

        return results_df

    def extract_hgmd_indel_coordinates(self):
        hgvs_hgmd_indel_1, hgvs_hgmd_indel_2, hgvs_hgmd_indel_3 = self.get_hgmd_indel_hgvs_dataframes()

        results_df = pd.concat([self[self.disease_column_name],                 # 1 column
                                self[self.gene_column_name],                    # 1 column
                                hgvs_hgmd_indel_3[2],                           # 1 column
                                hgvs_hgmd_indel_1[0],                           # 1 column
                                hgvs_hgmd_indel_1[1],                           # 1 column
                                hgvs_hgmd_indel_1[3],                           # 1 column
                                hgvs_hgmd_indel_2[3],                           # 1 column
                                hgvs_hgmd_indel_1[5],                           # 1 column
                                hgvs_hgmd_indel_2[5],                           # 1 column
                                self.get_changed_codon_position_dataframe(),    # 1 column
                                ],
                               axis=1
                               )

        results_df.columns = ["disease",
                              "gene",
                              "indel_type",
                              "indel_start",
                              "indel_end",
                              "del_sequence",
                              "del_length",
                              "ins_sequence",
                              "ins_length",
                              "codon_pos"]

        return results_df

    def extract_variant_coordinates(self, variant_type):
        if variant_type == "snp":
            return self.extract_snp_coordinates()
        elif variant_type == "indel":
            return self.extract_indel_coordinates()
        elif variant_type == "hgmd_indel":
            return self.extract_hgmd_indel_coordinates()
        else:
            raise ValueError("ERROR!!! Unknown or not implemented for coordinate extraction variant type")