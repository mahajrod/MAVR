#!/usr/bin/env python
import os
import time

import numpy as np

from Bio import SeqIO, Entrez

from Routines import FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict

from Tools.Abstract import Tool

class NCBIRoutines:
    def __init__(self):

        pass

    @staticmethod
    def efetch(database, id_list, out_file, retmode=None, rettype=None, seq_start=None, seq_stop=None, strand=None, verbose=False):
        # replacement for Biopython Entrez.efetch
        # Biopython Entrez.efetch is bugged - it ignores seq_start and seq_stop values
        # eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=669632474&retmode=text&rettype=gb&seq_start=10832751&seq_stop=10848091&strand=1
        query = "eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        query += "db=%s" % database
        query += "&id=%s" % (id_list if isinstance(id_list, str) else ",".join(id_list))
        query += "&retmode=%s" % retmode if retmode else ""
        query += "&rettype=%s" % rettype if rettype else ""

        query += "&seq_start=%s" % str(seq_start) if seq_start else ""
        query += "&seq_stop=%s" % str(seq_stop) if seq_stop else ""
        query += "&strand=%s" % str(strand) if strand else ""

        curl_string = "curl '%s' > %s" % (query, out_file)
        if verbose:
            print curl_string
        os.system(curl_string)

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

    def get_cds_for_proteins(self, protein_id_list, output_prefix, download_chunk_size=100, temp_dir="temp"):

        number_of_ids = len(protein_id_list)
        print "Totaly %i ids" % number_of_ids
        FileRoutines.save_mkdir(temp_dir)
        pep_file = "%s.pep.genbank" % output_prefix
        transcript_file = "%s.trascript.genbank" % output_prefix

        ranges = np.append(np.arange(0, number_of_ids, download_chunk_size), [number_of_ids])

        for i in range(0, len(ranges)-1):
            print "Downloading chunk %i" % i
            pep_tmp_file = "%s/%s_%i" % (temp_dir, pep_file, i)
            self.efetch("protein", protein_id_list[ranges[i]:ranges[i+1]], pep_tmp_file, rettype="gb", retmode="text")

        os.system("cat %s/* > %s" % (temp_dir, pep_file))
        print "BBBB"
        peptide_dict = SeqIO.index_db("tmp.idx", pep_file, format="genbank")
        downloaded_protein_ids = IdList(peptide_dict.keys())

        print "%i proteins were downloaded" % len(downloaded_protein_ids)
        not_downloded_proteins_ids = Tool.intersect_ids(protein_id_list, downloaded_protein_ids, mode="only_a")
        print "%i proteins were not downloaded" % len(not_downloded_proteins_ids)
        not_downloded_proteins_ids.write("%s.not_downloaded.ids" % output_prefix)


        pep_to_transcript_accordance = SynDict()
        for pep_id in peptide_dict:
            for feature in peptide_dict[pep_id]:
                if feature.type == "CDS":
                    print feature

    def get_cds_for_proteins_from_id_file(self, protein_id_file, output_prefix):
        pep_ids = IdList()
        pep_ids.read(protein_id_file)

        self.get_cds_for_proteins(pep_ids, output_prefix)
