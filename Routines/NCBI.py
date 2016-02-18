#!/usr/bin/env python
import os
import time

from Bio import SeqIO, Entrez

from Routines import FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict


class NCBIRoutines:
    def __init__(self):

        pass

    @staticmethod
    def efetch(database, id_list, out_file, retmode=None, rettype=None, seq_start=None, seq_stop=None, strand=None):
        # replacement for Biopython Entrez.efetch
        # Biopython Entrez.efetch is bugged - it ignores seq_start and seq_stop values
        # eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1
        query = "eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        query += "db=%s" % database
        query += "&id=%s" % (id_list if isinstance(id_list, str) else ",".join(id_list))
        query += "&retmode=%s" % retmode if retmode else None
        query += "&rettype=%s" % rettype if rettype else None

        query += "&seq_start=%s" % str(seq_start) if seq_start else None
        query += "&seq_stop=%s" % str(seq_stop) if seq_stop else None
        query += "&strand=%s" % str(strand) if strand else None

        os.system("curl '%s' > %s" % (query, out_file))

    def get_gene_sequences(self, email, query, retmax=100000, output_directory=None):
        if output_directory:
            FileRoutines.save_mkdir(output_directory)
        Entrez.email = email
        handle = Entrez.esearch(term=query, db="gene", rettype="xml", retmax=retmax)
        records = Entrez.read(handle, validate=False)
        gene_ids = records['IdList']
        for gene_id in gene_ids:
            print("Extracting sequence for gi %s" % gene_id)
            id_handle = Entrez.efetch("gene", id=gene_id, rettype="xml", retmax=retmax)
            record = Entrez.read(id_handle, validate=False)

            #print(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"])
            if "Entrezgene_locus" not in record[0]:
                continue
            if "Gene-commentary_seqs" not in record[0]["Entrezgene_locus"][0]:
                continue
            if "Seq-loc_int" not in record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]:
                continue
            gi = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_id']['Seq-id']['Seq-id_gi']
            start = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_from'])
            end = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_to']) + 1
            strand = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_strand']['Na-strand'].attributes['value']
            strand = 1 if strand == "plus" else 2
            # print(gi)
            # print(start, end, strand)
            # Biopython Entrez.efetch is bugged - it ignores seq_start and seq_stop values
            # eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1

            out_file = "%s/%s.gb" % (output_directory, gi)

            self.efetch("nuccore", gi, out_file, retmode="text", rettype="gb", seq_start=start, seq_stop=end,
                        strand=strand)
            time.sleep(0.4)
