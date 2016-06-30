#!/usr/bin/env python
import os
import time

import numpy as np

from Bio import SeqIO, Entrez

from Routines import FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict


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

    def get_cds_for_proteins(self, protein_id_list, output_prefix, download_chunk_size=100, temp_dir_prefix="temp"):

        from Tools.Abstract import Tool

        transcript_temp_dir = "%s_transcripts" % temp_dir_prefix
        protein_temp_dir = "%s_proteins" % temp_dir_prefix
        number_of_ids = len(protein_id_list)
        print "Total %i ids" % number_of_ids

        for directory in transcript_temp_dir, protein_temp_dir:
            FileRoutines.save_mkdir(directory)
        pep_file = "%s.pep.genbank" % output_prefix
        transcript_file = "%s.trascript.genbank" % output_prefix

        ranges = np.append(np.arange(0, number_of_ids, download_chunk_size), [number_of_ids])

        """
        print "Downloading proteins..."
        for i in range(0, len(ranges)-1):
            print "Downloading chunk %i" % i
            pep_tmp_file = "%s/%s_%i" % (protein_temp_dir, pep_file, i)
            self.efetch("protein", protein_id_list[ranges[i]:ranges[i+1]], pep_tmp_file, rettype="gb", retmode="text")

        os.system("cat %s/* > %s" % (protein_temp_dir, pep_file))
        """
        peptide_dict = SeqIO.index_db("tmp.idx", pep_file, format="genbank")
        downloaded_protein_ids = IdList(peptide_dict.keys())

        print "%i proteins were downloaded" % len(downloaded_protein_ids)
        not_downloded_proteins_ids = Tool.intersect_ids(protein_id_list, downloaded_protein_ids, mode="only_a")
        print "%i proteins were not downloaded" % len(not_downloded_proteins_ids)
        not_downloded_proteins_ids.write("%s.not_downloaded.ids" % output_prefix)

        pep_without_transcripts = IdList()
        pep_with_several_CDS_features = IdList()
        pep_to_transcript_accordance = SynDict()
        transcript_ids = IdList()

        print "Extracting transcript ids corresponding to proteins..."
        for pep_id in peptide_dict:
            for feature in peptide_dict[pep_id].features:
                if feature.type == "CDS":
                    try:
                        transcript_id = feature.qualifiers["coded_by"][0].split(":")[0]
                        if pep_id not in pep_to_transcript_accordance:
                            pep_to_transcript_accordance[pep_id] = [transcript_id]
                        else:


                            pep_to_transcript_accordance[pep_id].append(transcript_id)
                            print("Genbank record for %s contains several CDS features" % pep_id)
                            pep_with_several_CDS_features.append(pep_id)
                        if transcript_id in transcript_ids:
                            print "Repeated transcript id: %s" % transcript_id
                            continue
                        transcript_ids.append(transcript_id)
                    except:
                        print "Transcript id for %s was not found" % pep_id
                        pep_without_transcripts.append(pep_id)

        pep_with_several_CDS_features.write("%s.pep_with_several_CDS.ids" % output_prefix)
        pep_without_transcripts.write("%s.pep_without_transcripts.ids" % output_prefix)
        transcript_ids.write("%s.transcripts.ids" % output_prefix)

        number_of_transcripts = len(transcript_ids)
        print "%i transcripts were found" % number_of_transcripts

        pep_to_transcript_accordance.write("%s.pep_to_transcript.accordance" % output_prefix)

        transcript_ranges = np.append(np.arange(0, number_of_transcripts, download_chunk_size), [number_of_transcripts])

        print "Downloading transcripts..."
        for i in range(0, len(transcript_ranges)-1):
            print "Downloading chunk %i" % i
            transcript_tmp_file = "%s/%s_%i" % (transcript_temp_dir, transcript_file, i)
            self.efetch("nuccore", transcript_ids[transcript_ranges[i]:transcript_ranges[i+1]],
                        transcript_tmp_file, rettype="gb", retmode="text")

        os.system("cat %s/* > %s" % (transcript_temp_dir, transcript_file))


        transcript_dict = SeqIO.index_db("tmp_1.idx", transcript_file, format="genbank")
        """
        for transcript_id in transcript_dict:
            for feature in transcript_dict[transcript_id]:
                if feature.type == "CDS":
                    print feature
        """

        """
        for filename in "tmp.idx", "tmp_2.idx":
            os.remove(filename)

        """

    def get_cds_for_proteins_from_id_file(self, protein_id_file, output_prefix):
        pep_ids = IdList()
        pep_ids.read(protein_id_file)

        self.get_cds_for_proteins(pep_ids, output_prefix)
