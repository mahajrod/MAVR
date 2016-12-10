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
from Routines import FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict, TwoLvlDict


class AssemblySummary(OrderedDict):

    def __init__(self, entrez_summary_biopython):
        OrderedDict.__init__(self)
        parameters_dict = entrez_summary_biopython['DocumentSummarySet']['DocumentSummary'][0]
        self.stats = self.parse_stats_from_assembly_summary_metadata(parameters_dict['Meta'])

        for parameter in 'SpeciesName', 'Organism', 'Taxid', 'AssemblyName', 'FtpPath_GenBank', 'PropertyList', \
                         'SubmissionDate', 'AssemblyAccession', 'LastMajorReleaseAccession', \
                         'AssemblyType', 'PartialGenomeRepresentation', 'AssemblyClass', 'AnomalousList',\
                         'AssemblyStatus', 'BioSampleId':
            if parameter in parameters_dict:
                setattr(self, parameter, parameters_dict[parameter])
            else:
                setattr(self, parameter, None)
                #exec("self[%s] = None" % parameter)
        for entry in self.stats:
            self[entry] = self.stats[entry]
        for entry in parameters_dict:
            self[entry] = parameters_dict[entry]

        for parameter in "alt_loci_count", "chromosome_count", "contig_count", "contig_l50", "contig_n50", \
                         "non_chromosome_replicon_count", "replicon_count", "scaffold_count_all", \
                         "scaffold_count_placed", "scaffold_count_unlocalized", "scaffold_count_unplaced",\
                         "scaffold_l50", "scaffold_n50", "total_length", "ungapped_length":
            setattr(self, parameter, self[parameter])

        self.complete_genome = True if self.AssemblyStatus == 'Complete Genome' else False
        self.chromosome_lvl = True if 'has-chromosome' in self.PropertyList else False
        self.chloroplast = True if 'has-chloroplast' in self.PropertyList else False
        self.mitochondrion = True if 'has-mitochondrion' in self.PropertyList else False
        self.plasmid = True if 'has-plasmid' in self.PropertyList else False

    @staticmethod
    def parse_stats_from_assembly_summary_metadata(metadata):

        end = re.search("/Stats>", metadata).end()
        tmp_dict = xmltodict.parse(metadata[:end])
        stat_dict = OrderedDict()
        for entry in tmp_dict["Stats"]["Stat"]:
            if "scaffold_count" in entry['@category']:
                stat_dict[entry['@category'] + "_" + entry['@sequence_tag']] = int(entry['#text'])
            else:
                stat_dict[entry['@category']] = int(entry['#text'])

        return stat_dict

    def string_form(self, separator="\t"):
        pass


class NCBIRoutines:
    def __init__(self):

        pass

    @staticmethod
    def efetch(database, id_list, out_file, retmode=None, rettype=None, seq_start=None, seq_stop=None, strand=None, verbose=False,
               number_of_retries=100, retry_delay=None, log_file="efetch.log"):
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

        curl_options = " --retry %i" % number_of_retries if number_of_retries else ""
        curl_options += " --retry-delay %i" % retry_delay if retry_delay else ""
        curl_options += " '%s'" % query
        curl_options += " > %s" % out_file

        curl_string = "curl %s" % curl_options
        with open(log_file, "a") as log_fd:
            log_fd.write(curl_string + "\n")
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

        print "Downloading proteins..."
        for i in range(0, len(ranges)-1):
            print "Downloading chunk %i" % i
            pep_tmp_file = "%s/%s_%i" % (protein_temp_dir, pep_file, i)
            self.efetch("protein", protein_id_list[ranges[i]:ranges[i+1]], pep_tmp_file, rettype="gb", retmode="text")

        os.system("cat %s/* > %s" % (protein_temp_dir, pep_file))

        peptide_dict = SeqIO.index_db("tmp.idx", pep_file, format="genbank")
        downloaded_protein_ids = IdList(peptide_dict.keys())

        print "%i proteins were downloaded" % len(downloaded_protein_ids)
        not_downloaded_proteins_ids = Tool.intersect_ids([protein_id_list], [downloaded_protein_ids], mode="only_a")
        print "%i proteins were not downloaded" % len(not_downloaded_proteins_ids)
        not_downloaded_proteins_ids.write("%s.not_downloaded.ids" % output_prefix)
        downloaded_protein_ids.write("%s.downloaded.ids" % output_prefix)
        print Tool.intersect_ids([protein_id_list], [downloaded_protein_ids], mode="count")

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

        pep_to_transcript_accordance.write("%s.pep_to_transcript.accordance" % output_prefix, splited_values=True)

        transcript_ranges = np.append(np.arange(0, number_of_transcripts, download_chunk_size), [number_of_transcripts])

        print "Downloading transcripts..."
        for i in range(0, len(transcript_ranges)-1):
            print "Downloading chunk %i" % i
            transcript_tmp_file = "%s/%s_%i" % (transcript_temp_dir, transcript_file, i)
            self.efetch("nuccore", transcript_ids[transcript_ranges[i]:transcript_ranges[i+1]],
                        transcript_tmp_file, rettype="gb", retmode="text")

        os.system("cat %s/* > %s" % (transcript_temp_dir, transcript_file))


        transcript_dict = SeqIO.index_db("tmp_1.idx", transcript_file, format="genbank")

        cds_records_list = []
        for transcript_id in transcript_dict:
            for feature in transcript_dict[transcript_id].features:
                CDS_counter = 1
                if feature.type == "CDS":
                    #print feature

                    feature_seq = feature.extract(transcript_dict[transcript_id].seq)
                    feature_id = transcript_id  # case with several CDS per transcripts is was not taken into account
                    if "protein_id" in feature.qualifiers:
                        description = "protein=%s" % feature.qualifiers["protein_id"][0]
                    else:
                        print "Corresponding protein id was not found for %s" % transcript_id
                    cds_records_list.append(SeqRecord(seq=feature_seq, id=feature_id, description=description))
        SeqIO.write(cds_records_list, "%s.cds" % output_prefix, format="fasta")

        stat_string = "Input protein ids\t %i\n" % number_of_ids
        stat_string += "Downloaded proteins\t%i\n" % number_of_transcripts
        stat_string += "Downloaded transcripts\t%i\n" % len(transcript_dict)

        print stat_string

        with open("%s.stats" % output_prefix, "w") as stat_fd:
            stat_fd.write(stat_string)

        for filename in "tmp.idx", "tmp_1.idx":
            os.remove(filename)

    def get_cds_for_proteins_from_id_file(self, protein_id_file, output_prefix):
        pep_ids = IdList()
        pep_ids.read(protein_id_file)

        self.get_cds_for_proteins(pep_ids, output_prefix)

    @staticmethod
    def get_taxonomy(taxa_list, output_file, email, input_type="latin"):
        Entrez.email = email
        out_file = open(output_file, "w")
        out_file.write("#species\tlineage\n")
        if input_type == "latin":
            for taxon in taxa_list:
                print "Handling %s" % taxon
                summary = Entrez.read(Entrez.esearch(db="taxonomy", term=taxon))
                if summary:
                    id_list = summary["IdList"]
                    for id in id_list:
                        print "handling %s" % id
                        record = Entrez.read(Entrez.efetch(db="taxonomy", id=id, retmode="xml"))
                        out_file.write("%s\t%s\t%s\n" % (taxon, record[0]["Rank"], record[0]["Lineage"]))

        elif input_type == "id":
            for taxon in taxa_list:
                print "Handling %s" % taxon
                record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxon, retmode="xml"))
                out_file.write("%s\t%s\t%s\n" % (taxon, record[0]["Rank"], record[0]["Lineage"]))

    def get_taxonomy_from_id_file(self, taxa_file, output_file, email, input_type="latin"):

        taxa_list = IdList()
        taxa_list.read(taxa_file)

        self.get_taxonomy(taxa_list, output_file, email, input_type=input_type)

    def get_taxa_genomes_summary(self, taxa, email, output_file):
        Entrez.email = email
        taxa_list = taxa if isinstance(taxa, Iterable) else [taxa]

        for taxon in taxa_list:
            search_term = "%s[Orgn]" % taxon

            summary = Entrez.read(Entrez.esearch(db="genome", term=search_term, retmax=10000, retmode="xml"))
            print "Were found %s species" % summary["Count"]
            #print summary

            for species_id in summary["IdList"]: #[167] :
                print "Species id %s " % species_id

                species_summary = Entrez.read(Entrez.esummary(db="genome", id=species_id, retmax=10000, retmode="xml"))
                #print species_summary

                assembly_links = Entrez.read(Entrez.elink(dbfrom="genome", id=species_id, retmode="xml",
                                                          retmax=10000, linkname="genome_assembly"))
                #print len(links)
                #print links
                #print links[0]["LinkSetDb"][0]["Link"]
                assembly_ids = [id_dict["Id"] for id_dict in assembly_links[0]["LinkSetDb"][0]["Link"]]
                #print len(assembly_links[0]["LinkSetDb"][0]["Link"])
                #print assembly_ids
                #print len(assembly_ids)
                assembly_dict = TwoLvlDict()
                assemblies_with_ambiguous_taxonomies = SynDict()
                for assembly_id in assembly_ids: #[:1]:
                    assembly_summary = AssemblySummary(Entrez.read(Entrez.esummary(db="assembly",
                                                                                   id=assembly_id,
                                                                                   retmode="xml")))
                    print assembly_summary
                    #print assembly_summary.__dict__

                    print assembly_summary['PropertyList']
                    print assembly_summary.SpeciesName
                    """
                                                    u'ChainId': '1886755',
                                                    u'AsmUpdateDate': '2016/11/28 00:00',
                                                    u'GbUid': '3745278',
                                                    u'PropertyList': ['full-genome-representation', 'has-chromosome', 'has-plasmid', 'latest', 'latest_genbank'],
                                                    u'SubmissionDate': '2016/11/28 00:00',
                                                    u'GB_BioProjects': [{u'BioprojectAccn': 'PRJNA353939',
                                                                         u'BioprojectId': '353939'}],
                                                    u'AssemblyAccession': 'GCA_001886755.1',
                                                    u'LastMajorReleaseAccession': 'GCA_001886755.1',
                                                    u'Synonym': {u'RefSeq': '',
                                                                 u'Genbank': 'GCA_001886755.1',
                                                                 u'Similarity': ''},
                                                    u'FtpPath_RefSeq': '',
                                                    u'AsmReleaseDate_RefSeq': '1/01/01 00:00',
                                                    u'AssemblyType': 'haploid',
                                                    u'AssemblyDescription': '',
                                                    u'AsmReleaseDate_GenBank': '2016/11/28 00:00',
                                                    u'RS_Projects': [],
                                                    u'RefSeq_category': 'na',
                                                    u'PartialGenomeRepresentation': 'false',
                                                    u'Coverage': '0',
                                                    u'AssemblyClass': 'haploid',
                                                    u'ExclFromRefSeq': [],
                                                    u'SeqReleaseDate': '2016/11/28 00:00',
                                                    u'AnomalousList': [],
                                                    u'AssemblyStatus': 'Complete Genome',
                                                    u'AsmReleaseDate': '2016/11/28 00:00',
                                                    u'ReleaseLevel': 'Major',
                                                    u'Taxid': '562',
                                                    u'LastUpdateDate': '2016/11/28 00:00',
                                                    u'RsUid': '',
                                                    u'FromType': '',
                                                    u'WGS': '',
                                                    u'GB_Projects': ['353939'],
                                                    u'BioSampleId': '6029894',
                                                    u'AssemblyName': 'ASM188675v1',
                                                    u'EnsemblName': '',
                                                    u'Organism': 'Escherichia coli (E. coli)',
                                                    u'Biosource': {u'Isolate': '',
                                                                   u'InfraspeciesList': [{u'Sub_type': 'strain',
                                                                                          u'Sub_value': 'MRSN346595'}],
                                                                   u'Sex': ''},
                                                    u'SpeciesName': 'Escherichia coli',
                                                    u'BioSampleAccn': 'SAMN06029894',
                                                    u'FtpPath_GenBank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1',
                                                    u'SpeciesTaxid': '562',
                                                    u'UCSCName': '',
                                                    u'Primary': '3745268',
                                                    u'Meta': ' <Stats> <Stat category="alt_loci_count" sequence_tag="all">0</Stat> <Stat category="chromosome_count" sequence_tag="all">1</Stat> <Stat category="contig_count" sequence_tag="all">6</Stat> <Stat category="contig_l50" sequence_tag="all">1</Stat> <Stat category="contig_n50" sequence_tag="all">4796421</Stat> <Stat category="non_chromosome_replicon_count" sequence_tag="all">5</Stat> <Stat category="replicon_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="placed">6</Stat> <Stat category="scaffold_count" sequence_tag="unlocalized">0</Stat> <Stat category="scaffold_count" sequence_tag="unplaced">0</Stat> <Stat category="scaffold_l50" sequence_tag="all">1</Stat> <Stat category="scaffold_n50" sequence_tag="all">4796421</Stat> <Stat category="total_length" sequence_tag="all">5116627</Stat> <Stat category="ungapped_length" sequence_tag="all">5116627</Stat> </Stats> <FtpSites>   <FtpPath type="GenBank">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1</FtpPath> </FtpSites> <assembly-level>5</assembly-level> <assembly-status>Complete Genome</assembly-status> <representative-status>na</representative-status> <submitter-organization>Walter Reed Army Institute of Research</submitter-organization>    ',
                                                    u'NCBIReleaseDate': '2016/11/28 00:00',
                                                    u'SubmitterOrganization': 'Walter Reed Army Institute of Research',
                                                    u'RS_BioProjects': [],
                                                    u'SortOrder': '5C10018867559899'},
                                                    attributes={u'uid': u'894571'}
                    """
                    """
                    taxonomy_links = Entrez.read(Entrez.elink(dbfrom="assembly", id=assembly_id, retmode="xml",
                                                              retmax=10000, linkname="assembly_taxonomy"))
                    #print taxonomy_links
                    taxonomy_ids = [id_dict["Id"] for id_dict in taxonomy_links[0]["LinkSetDb"][0]["Link"]]

                    if len(taxonomy_ids) > 1:
                        print "WARNING: more than one taxonomy id for assembly %s" % assembly_ids
                        assemblies_with_ambiguous_taxonomies[assembly_id] = taxonomy_ids
                    #print taxonomy_ids

                    taxonomy_record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxonomy_ids[0], retmode="xml"))

                    #print taxonomy_record
                    scientific_name = taxonomy_record[0]["ScientificName"]
                    lineage = taxonomy_record[0]["Lineage"]
                    #common_name_list = taxonomy_record[0]["CommonName"]
                    print taxonomy_record[0]["ScientificName"]
                    """
                print "\n\n"
        pass
"""
assembly_nuccore 	Nucleotide 	Nucleotide 	Nucleotide 	5000
assembly_nuccore_insdc 	Nucleotide INSDC 	Nucleotide INSDC 	Nucleotide INSDC (GenBank) 	5000
assembly_nuccore_refseq 	Nucleotide RefSeq 	Nucleotide RefSeq 	Nucleotide RefSeq 	5000
assembly_nuccore_wgsmaster 	WGS Master 	WGS Master 	WGS Master 	5000
assembly_nuccore_wgscontig 	WGS contigs 	WGS contigs 	WGS contigs 	5000
assembly_nuccore_insdc 	Nucleotide INSDC 	Nucleotide INSDC 	Nucleotide INSDC (GenBank) 	5000
assembly_nuccore_refseq 	Nucleotide RefSeq 	Nucleotide RefSeq 	Nucleotide RefSeq 	5000
assembly_nuccore_wgsmaster 	WGS Master 	WGS Master 	WGS Master
"""

"""

#Example of assembly summary

{u'DocumentSummarySet':
     DictElement({u'DbBuild': 'Build161130-1600.1',
                  u'DocumentSummary': [DictElement({u'ChainId': '1886755',
                                                    u'AsmUpdateDate': '2016/11/28 00:00',
                                                    u'GbUid': '3745278',
                                                    u'PropertyList': ['full-genome-representation', 'has-chromosome', 'has-plasmid', 'latest', 'latest_genbank'],
                                                    u'SubmissionDate': '2016/11/28 00:00',
                                                    u'GB_BioProjects': [{u'BioprojectAccn': 'PRJNA353939',
                                                                         u'BioprojectId': '353939'}],
                                                    u'AssemblyAccession': 'GCA_001886755.1',
                                                    u'LastMajorReleaseAccession': 'GCA_001886755.1',
                                                    u'Synonym': {u'RefSeq': '',
                                                                 u'Genbank': 'GCA_001886755.1',
                                                                 u'Similarity': ''},
                                                    u'FtpPath_RefSeq': '',
                                                    u'AsmReleaseDate_RefSeq': '1/01/01 00:00',
                                                    u'AssemblyType': 'haploid',
                                                    u'AssemblyDescription': '',
                                                    u'AsmReleaseDate_GenBank': '2016/11/28 00:00',
                                                    u'RS_Projects': [],
                                                    u'RefSeq_category': 'na',
                                                    u'PartialGenomeRepresentation': 'false',
                                                    u'Coverage': '0',
                                                    u'AssemblyClass': 'haploid',
                                                    u'ExclFromRefSeq': [],
                                                    u'SeqReleaseDate': '2016/11/28 00:00',
                                                    u'AnomalousList': [],
                                                    u'AssemblyStatus': 'Complete Genome',
                                                    u'AsmReleaseDate': '2016/11/28 00:00',
                                                    u'ReleaseLevel': 'Major',
                                                    u'Taxid': '562',
                                                    u'LastUpdateDate': '2016/11/28 00:00',
                                                    u'RsUid': '',
                                                    u'FromType': '',
                                                    u'WGS': '',
                                                    u'GB_Projects': ['353939'],
                                                    u'BioSampleId': '6029894',
                                                    u'AssemblyName': 'ASM188675v1',
                                                    u'EnsemblName': '',
                                                    u'Organism': 'Escherichia coli (E. coli)',
                                                    u'Biosource': {u'Isolate': '',
                                                                   u'InfraspeciesList': [{u'Sub_type': 'strain',
                                                                                          u'Sub_value': 'MRSN346595'}],
                                                                   u'Sex': ''},
                                                    u'SpeciesName': 'Escherichia coli',
                                                    u'BioSampleAccn': 'SAMN06029894',
                                                    u'FtpPath_GenBank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1',
                                                    u'SpeciesTaxid': '562',
                                                    u'UCSCName': '',
                                                    u'Primary': '3745268',
                                                    u'Meta': ' <Stats> <Stat category="alt_loci_count" sequence_tag="all">0</Stat> <Stat category="chromosome_count" sequence_tag="all">1</Stat> <Stat category="contig_count" sequence_tag="all">6</Stat> <Stat category="contig_l50" sequence_tag="all">1</Stat> <Stat category="contig_n50" sequence_tag="all">4796421</Stat> <Stat category="non_chromosome_replicon_count" sequence_tag="all">5</Stat> <Stat category="replicon_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="all">6</Stat> <Stat category="scaffold_count" sequence_tag="placed">6</Stat> <Stat category="scaffold_count" sequence_tag="unlocalized">0</Stat> <Stat category="scaffold_count" sequence_tag="unplaced">0</Stat> <Stat category="scaffold_l50" sequence_tag="all">1</Stat> <Stat category="scaffold_n50" sequence_tag="all">4796421</Stat> <Stat category="total_length" sequence_tag="all">5116627</Stat> <Stat category="ungapped_length" sequence_tag="all">5116627</Stat> </Stats> <FtpSites>   <FtpPath type="GenBank">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/755/GCA_001886755.1_ASM188675v1</FtpPath> </FtpSites> <assembly-level>5</assembly-level> <assembly-status>Complete Genome</assembly-status> <representative-status>na</representative-status> <submitter-organization>Walter Reed Army Institute of Research</submitter-organization>    ',
                                                    u'NCBIReleaseDate': '2016/11/28 00:00',
                                                    u'SubmitterOrganization': 'Walter Reed Army Institute of Research',
                                                    u'RS_BioProjects': [],
                                                    u'SortOrder': '5C10018867559899'},
                                                    attributes={u'uid': u'894571'}
                                                   )
                                       ]
                  },
                 attributes={u'status': u'OK'}
                 )
 }
"""
