#!/usr/bin/env python
import os

from Bio import SeqIO, Entrez


from CustomCollections.GeneralCollections import IdList, SynDict


class NCBIRoutines:
    def __init__(self):

        pass

    @staticmethod
    def get_gene_sequences(email, query, retmax=100000):
        Entrez.email = email
        handle = Entrez.esearch(term=query, db="gene", rettype="xml", retmax=retmax)
        records = Entrez.read(handle, validate=False)
        gene_ids = records['IdList']
        for gene_id in gene_ids[0:1]:
            print(gene_id)
            id_handle = Entrez.efetch("gene", id=gene_id, rettype="xml", retmax=retmax)
            record = Entrez.read(id_handle, validate=False)

            print(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"])

            gi = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_id']['Seq-id']['Seq-id_gi']
            start = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_from'])
            end = int(record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_to'])
            strand = record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]['Seq-interval_strand']['Na-strand'].attributes['value']
            strand = 1 if strand == "plus" else -1
            print(gi)
            print(start, end, strand)
            #eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1
            gene_handle = Entrez.efetch(db="nuccore",
                                        id=gi,
                                        rettype="gb",
                                        retmode="text",
                                        seq_start=start,
                                        seq_stop=strand,
                                        strand=1)
            record = SeqIO.read(gene_handle, "genbank")
            print(record)
            SeqIO.write(record, "tmp.gb", format="genbank")
            #print(entry["Entrezgene_locus"][0].keys())
            #print record.keys()
