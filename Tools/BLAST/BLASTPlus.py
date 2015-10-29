#!/usr/bin/env python
__author__ = 'mahajrod'

import os
from subprocess import PIPE, Popen

from Tools.Abstract import Tool
from Tools.LinuxTools import CGAS
from Routines.File import check_path, save_mkdir, split_filename


class BLASTPlus(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blastn", path=path, max_threads=max_threads)

    def parallel_blast(self, blast_command, seqfile, database, outfile=None,
                       blast_options=None, split_dir="splited_fasta",
                       splited_output_dir="splited_output_dir",
                       evalue=None, output_format=None,
                       threads=None, num_of_seqs_per_scan=None,
                       combine_output_to_single_file=True):

        splited_dir = check_path(split_dir)
        splited_out_dir = check_path(splited_output_dir)
        save_mkdir(splited_dir)
        save_mkdir(splited_out_dir)

        number_of_files = num_of_seqs_per_scan if num_of_seqs_per_scan else 5 * threads if threads else 5 * self.threads
        self.split_fasta(seqfile, splited_dir, num_of_files=number_of_files)
        input_list_of_files = sorted(os.listdir(splited_dir))
        list_of_files = []

        for filename in input_list_of_files:
            filename_prefix = split_filename(filename)[1]

            input_file = "%s%s" % (splited_dir, filename)
            output_file = "%s%s.hits" % (splited_out_dir, filename_prefix)

            list_of_files.append((input_file, output_file))

        options_list = []
        out_files = []

        for in_file, out_filename in list_of_files:

            options = " -out %s" % out_filename

            options += " -db %s" % database
            options += " -query %s" % seqfile
            options += " %s" % blast_options
            options += " -evalue %f" % evalue if evalue else ""
            options += " -outfmt %i" % output_format if output_format else ""
            options_list.append(options)
            out_files.append(out_filename)

        self.parallel_execute(options_list, cmd=blast_command, threads=threads)

        if combine_output_to_single_file:
            CGAS.cat(out_files, output=outfile)


class BLASTn(Tool, BLASTPlus):
    """
    USAGE
      blastn [-h] [-help] [-import_search_strategy filename]
        [-export_search_strategy filename] [-task task_name] [-db database_name]
        [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
        [-negative_gilist filename] [-entrez_query entrez_query]
        [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
        [-subject subject_input_file] [-subject_loc range] [-query input_file]
        [-out output_file] [-evalue evalue] [-word_size int_value]
        [-gapopen open_penalty] [-gapextend extend_penalty]
        [-perc_identity float_value] [-xdrop_ungap float_value]
        [-xdrop_gap float_value] [-xdrop_gap_final float_value]
        [-searchsp int_value] [-max_hsps_per_subject int_value] [-penalty penalty]
        [-reward reward] [-no_greedy] [-min_raw_gapped_score int_value]
        [-template_type type] [-template_length int_value] [-dust DUST_options]
        [-filtering_db filtering_database]
        [-window_masker_taxid window_masker_taxid]
        [-window_masker_db window_masker_db] [-soft_masking soft_masking]
        [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
        [-best_hit_score_edge float_value] [-window_size int_value]
        [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
        [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
        [-outfmt format] [-show_gis] [-num_descriptions int_value]
        [-num_alignments int_value] [-html] [-max_target_seqs num_sequences]
        [-num_threads int_value] [-remote] [-version]

    DESCRIPTION
       Nucleotide-Nucleotide BLAST 2.2.28+

    OPTIONAL ARGUMENTS
     -h
       Print USAGE and DESCRIPTION;  ignore all other parameters
     -help
       Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
     -version
       Print version number;  ignore other arguments

     *** Input query options
     -query <File_In>
       Input file name
       Default = `-'
     -query_loc <String>
       Location on the query sequence in 1-based offsets (Format: start-stop)
     -strand <String, `both', `minus', `plus'>
       Query strand(s) to search against database/subject
       Default = `both'

     *** General search options
     -task <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast'
                    'megablast' 'rmblastn' >
       Task to execute
       Default = `megablast'
     -db <String>
       BLAST database name
        * Incompatible with:  subject, subject_loc
     -out <File_Out>
       Output file name
       Default = `-'
     -evalue <Real>
       Expectation value (E) threshold for saving hits
       Default = `10'
     -word_size <Integer, >=4>
       Word size for wordfinder algorithm (length of best perfect match)
     -gapopen <Integer>
       Cost to open a gap
     -gapextend <Integer>
       Cost to extend a gap
     -penalty <Integer, <=0>
       Penalty for a nucleotide mismatch
     -reward <Integer, >=0>
       Reward for a nucleotide match
     -use_index <Boolean>
       Use MegaBLAST database index
     -index_name <String>
       MegaBLAST database index name

     *** BLAST-2-Sequences options
     -subject <File_In>
       Subject sequence(s) to search
        * Incompatible with:  db, gilist, seqidlist, negative_gilist,
       db_soft_mask, db_hard_mask
     -subject_loc <String>
       Location on the subject sequence in 1-based offsets (Format: start-stop)
        * Incompatible with:  db, gilist, seqidlist, negative_gilist,
       db_soft_mask, db_hard_mask, remote

     *** Formatting options
     -outfmt <String>
       alignment view options:
         0 = pairwise,
         1 = query-anchored showing identities,
         2 = query-anchored no identities,
         3 = flat query-anchored, show identities,
         4 = flat query-anchored, no identities,
         5 = XML Blast output,
         6 = tabular,
         7 = tabular with comment lines,
         8 = Text ASN.1,
         9 = Binary ASN.1,
        10 = Comma-separated values,
        11 = BLAST archive format (ASN.1)

       Options 6, 7, and 10 can be additionally configured to produce
       a custom format specified by space delimited format specifiers.
       The supported format specifiers are:
            qseqid means Query Seq-id
               qgi means Query GI
              qacc means Query accesion
           qaccver means Query accesion.version
              qlen means Query sequence length
            sseqid means Subject Seq-id
         sallseqid means All subject Seq-id(s), separated by a ';'
               sgi means Subject GI
            sallgi means All subject GIs
              sacc means Subject accession
           saccver means Subject accession.version
           sallacc means All subject accessions
              slen means Subject sequence length
            qstart means Start of alignment in query
              qend means End of alignment in query
            sstart means Start of alignment in subject
              send means End of alignment in subject
              qseq means Aligned part of query sequence
              sseq means Aligned part of subject sequence
            evalue means Expect value
          bitscore means Bit score
             score means Raw score
            length means Alignment length
            pident means Percentage of identical matches
            nident means Number of identical matches
          mismatch means Number of mismatches
          positive means Number of positive-scoring matches
           gapopen means Number of gap openings
              gaps means Total number of gaps
              ppos means Percentage of positive-scoring matches
            frames means Query and subject frames separated by a '/'
            qframe means Query frame
            sframe means Subject frame
              btop means Blast traceback operations (BTOP)
           staxids means Subject Taxonomy ID(s), separated by a ';'
         sscinames means Subject Scientific Name(s), separated by a ';'
         scomnames means Subject Common Name(s), separated by a ';'
        sblastnames means Subject Blast Name(s), separated by a ';'
                 (in alphabetical order)
        sskingdoms means Subject Super Kingdom(s), separated by a ';'
                 (in alphabetical order)
            stitle means Subject Title
        salltitles means All Subject Title(s), separated by a '<>'
           sstrand means Subject Strand
             qcovs means Query Coverage Per Subject
           qcovhsp means Query Coverage Per HSP
       When not provided, the default value is:
       'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
       evalue bitscore', which is equivalent to the keyword 'std'
       Default = `0'
     -show_gis
       Show NCBI GIs in deflines?
     -num_descriptions <Integer, >=0>
       Number of database sequences to show one-line descriptions for
       Not applicable for outfmt > 4
       Default = `500'
        * Incompatible with:  max_target_seqs
     -num_alignments <Integer, >=0>
       Number of database sequences to show alignments for
       Default = `250'
        * Incompatible with:  max_target_seqs
     -html
       Produce HTML output?

     *** Query filtering options
     -dust <String>
       Filter query sequence with DUST (Format: 'yes', 'level window linker', or
       'no' to disable)
       Default = `20 64 1'
     -filtering_db <String>
       BLAST database containing filtering elements (i.e.: repeats)
     -window_masker_taxid <Integer>
       Enable WindowMasker filtering using a Taxonomic ID
     -window_masker_db <String>
       Enable WindowMasker filtering using this repeats database.
     -soft_masking <Boolean>
       Apply filtering locations as soft masks
       Default = `true'
     -lcase_masking
       Use lower case filtering in query and subject sequence(s)?

     *** Restrict search or results
     -gilist <String>
       Restrict search of database to list of GI's
        * Incompatible with:  negative_gilist, seqidlist, remote, subject,
       subject_loc
     -seqidlist <String>
       Restrict search of database to list of SeqId's
        * Incompatible with:  gilist, negative_gilist, remote, subject,
       subject_loc
     -negative_gilist <String>
       Restrict search of database to everything except the listed GIs
        * Incompatible with:  gilist, seqidlist, remote, subject, subject_loc
     -entrez_query <String>
       Restrict search with the given Entrez query
        * Requires:  remote
     -db_soft_mask <String>
       Filtering algorithm ID to apply to the BLAST database as soft masking
        * Incompatible with:  db_hard_mask, subject, subject_loc
     -db_hard_mask <String>
       Filtering algorithm ID to apply to the BLAST database as hard masking
        * Incompatible with:  db_soft_mask, subject, subject_loc
     -perc_identity <Real, 0..100>
       Percent identity
     -culling_limit <Integer, >=0>
       If the query range of a hit is enveloped by that of at least this many
       higher-scoring hits, delete the hit
        * Incompatible with:  best_hit_overhang, best_hit_score_edge
     -best_hit_overhang <Real, (>=0 and =<0.5)>
       Best Hit algorithm overhang value (recommended value: 0.1)
        * Incompatible with:  culling_limit
     -best_hit_score_edge <Real, (>=0 and =<0.5)>
       Best Hit algorithm score edge value (recommended value: 0.1)
        * Incompatible with:  culling_limit
     -max_target_seqs <Integer, >=1>
       Maximum number of aligned sequences to keep
       Not applicable for outfmt <= 4
       Default = `500'
        * Incompatible with:  num_descriptions, num_alignments

     *** Discontiguous MegaBLAST options
     -template_type <String, `coding', `coding_and_optimal', `optimal'>
       Discontiguous MegaBLAST template type
        * Requires:  template_length
     -template_length <Integer, Permissible values: '16' '18' '21' >
       Discontiguous MegaBLAST template length
        * Requires:  template_type

     *** Statistical options
     -dbsize <Int8>
       Effective length of the database
     -searchsp <Int8, >=0>
       Effective length of the search space
     -max_hsps_per_subject <Integer, >=0>
       Override maximum number of HSPs per subject to save for ungapped searches
       (0 means do not override)
       Default = `0'

     *** Search strategy options
     -import_search_strategy <File_In>
       Search strategy to use
        * Incompatible with:  export_search_strategy
     -export_search_strategy <File_Out>
       File name to record the search strategy used
        * Incompatible with:  import_search_strategy

     *** Extension options
     -xdrop_ungap <Real>
       X-dropoff value (in bits) for ungapped extensions
     -xdrop_gap <Real>
       X-dropoff value (in bits) for preliminary gapped extensions
     -xdrop_gap_final <Real>
       X-dropoff value (in bits) for final gapped alignment
     -no_greedy
       Use non-greedy dynamic programming extension
     -min_raw_gapped_score <Integer>
       Minimum raw gapped score to keep an alignment in the preliminary gapped and
       traceback stages
     -ungapped
       Perform ungapped alignment only?
     -window_size <Integer, >=0>
       Multiple hits window size, use 0 to specify 1-hit algorithm
     -off_diagonal_range <Integer, >=0>
       Number of off-diagonals to search for the 2nd hit, use 0 to turn off
       Default = `0'

     *** Miscellaneous options
     -parse_deflines
       Should the query and subject defline(s) be parsed?
     -num_threads <Integer, >=1>
       Number of threads (CPUs) to use in the BLAST search
       Default = `1'
        * Incompatible with:  remote
     -remote
       Execute search remotely?
        * Incompatible with:  gilist, seqidlist, negative_gilist, subject_loc,
       num_threads
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blastn", path=path, max_threads=max_threads)

    def parallel_blastn(self, seqfile, database, outfile=None,
                        blast_options=None, split_dir="splited_fasta",
                        splited_output_dir="splited_output_dir",
                        evalue=None, output_format=None,
                        threads=None, num_of_seqs_per_scan=None,
                        combine_output_to_single_file=True):

        self.parallel_blast("blastn", seqfile, database, outfile=outfile,
                            blast_options=blast_options, split_dir=split_dir,
                            splited_output_dir=splited_output_dir,
                            evalue=evalue, output_format=output_format,
                            threads=threads, num_of_seqs_per_scan=num_of_seqs_per_scan,
                            combine_output_to_single_file=combine_output_to_single_file)


class BLASTp(Tool, BLASTPlus):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blastp", path=path, max_threads=max_threads)

    def parallel_blastp(self, seqfile, database, outfile=None,
                        blast_options=None, split_dir="splited_fasta",
                        splited_output_dir="splited_output_dir",
                        evalue=None, output_format=None,
                        threads=None, num_of_seqs_per_scan=None,
                        combine_output_to_single_file=True):

        self.parallel_blast("blastp", seqfile, database, outfile=outfile,
                            blast_options=blast_options, split_dir=split_dir,
                            splited_output_dir=splited_output_dir,
                            evalue=evalue, output_format=output_format,
                            threads=threads, num_of_seqs_per_scan=num_of_seqs_per_scan,
                            combine_output_to_single_file=combine_output_to_single_file)


class DustMasker(Tool):
    def __init__(self, path="", max_threads=4):
        """
        USAGE
          dustmasker [-h] [-help] [-xmlhelp] [-in input_file_name]
            [-out output_file_name] [-window window_size] [-level level]
            [-linker linker] [-infmt input_format] [-outfmt output_format]
            [-parse_seqids] [-version-full]

        DESCRIPTION
           Low complexity region masker based on Symmetric DUST algorithm

        OPTIONAL ARGUMENTS
         -h
           Print USAGE and DESCRIPTION;  ignore all other parameters
         -help
           Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
         -xmlhelp
           Print USAGE, DESCRIPTION and ARGUMENTS in XML format; ignore all other
           parameters
         -in <File_In>
           input file name
           Default = `-'
         -out <File_Out>
           output file name
           Default = `-'
         -window <Integer>
           DUST window length
           Default = `64'
         -level <Integer>
           DUST level (score threshold for subwindows)
           Default = `20'
         -linker <Integer>
           DUST linker (how close masked intervals should be to get merged together).
           Default = `1'
         -infmt <String>
           input format (possible values: fasta, blastdb)
           Default = `fasta'
         -outfmt <String, `acclist', `fasta', `interval', `maskinfo_asn1_bin',
                          `maskinfo_asn1_text', `maskinfo_xml', `seqloc_asn1_bin',
                          `seqloc_asn1_text', `seqloc_xml'>
           output format
           Default = `interval'
         -parse_seqids
           Parse Seq-ids in FASTA input
         -version-full
           Print extended version data;  ignore other arguments
        """
        Tool.__init__(self, "dustmasker", path=path, max_threads=max_threads)

    def mask(self, input_file, output_file, output_format="maskinfo_asn1_bin", input_format="fasta", parse_seqids=True):
        options = " "
        options += " -in %s" % input_file
        options += " -infmt %s" % input_format
        options += " -outfmt %s" % output_format
        options += " -parse_seqids" if parse_seqids else ""
        options += " -out %s" % output_file

        self.execute(options)


class MakeBLASTDb(Tool):
    def __init__(self, path="", max_threads=4):
        """
        USAGE
          makeblastdb [-h] [-help] [-in input_file] [-input_type type]
            -dbtype molecule_type [-title database_title] [-parse_seqids]
            [-hash_index] [-mask_data mask_data_files] [-mask_id mask_algo_ids]
            [-mask_desc mask_algo_descriptions] [-gi_mask]
            [-gi_mask_name gi_based_mask_names] [-out database_name]
            [-max_file_sz number_of_bytes] [-logfile File_Name] [-taxid TaxID]
            [-taxid_map TaxIDMapFile] [-version]

        DESCRIPTION
           Application to create BLAST databases, version 2.2.30+

        REQUIRED ARGUMENTS
         -dbtype <String, `nucl', `prot'>
           Molecule type of target db

        OPTIONAL ARGUMENTS
         -h
           Print USAGE and DESCRIPTION;  ignore all other parameters
         -help
           Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
         -version
           Print version number;  ignore other arguments

         *** Input options
         -in <File_In>
           Input file/database name
           Default = `-'
         -input_type <String, `asn1_bin', `asn1_txt', `blastdb', `fasta'>
           Type of the data specified in input_file
           Default = `fasta'

         *** Configuration options
         -title <String>
           Title for BLAST database
           Default = input file name provided to -in argument
         -parse_seqids
           Option to parse seqid for FASTA input if set, for all other input types
           seqids are parsed automatically
         -hash_index
           Create index of sequence hash values.

         *** Sequence masking options
         -mask_data <String>
           Comma-separated list of input files containing masking data as produced by
           NCBI masking applications (e.g. dustmasker, segmasker, windowmasker)
         -mask_id <String>
           Comma-separated list of strings to uniquely identify the masking algorithm
            * Requires:  mask_data
            * Incompatible with:  gi_mask
         -mask_desc <String>
           Comma-separated list of free form strings to describe the masking algorithm
           details
            * Requires:  mask_id
         -gi_mask
           Create GI indexed masking data.
            * Requires:  parse_seqids
            * Incompatible with:  mask_id
         -gi_mask_name <String>
           Comma-separated list of masking data output files.
            * Requires:  mask_data, gi_mask

         *** Output options
         -out <String>
           Name of BLAST database to be created
           Default = input file name provided to -in argumentRequired if multiple
           file(s)/database(s) are provided as input
         -max_file_sz <String>
           Maximum file size for BLAST database files
           Default = `1GB'
         -logfile <File_Out>
           File to which the program log should be redirected

         *** Taxonomy options
         -taxid <Integer, >=0>
           Taxonomy ID to assign to all sequences
            * Incompatible with:  taxid_map
         -taxid_map <File_In>
           Text file mapping sequence IDs to taxonomy IDs.
           Format:<SequenceId> <TaxonomyId><newline>
            * Requires:  parse_seqids
            * Incompatible with:  taxid

        """
        Tool.__init__(self, "makeblastdb", path=path, max_threads=max_threads)

        #makes BLAST database from fasta file

    def make_db(self, input_file, db_title, mask_data, db_type, output_file=None, input_format="fasta", parse_seqids=True):
        # mask_data can be either string or list of strings
        options = " -dbtype %s" % db_type
        options += " -in %s" % input_file
        options += " -input_type %s" % input_format
        options += " -parse_seqids" if parse_seqids else ""
        options += " -title %s" % db_title
        options += " -out %s" % output_file if None else " -out %s" % db_title
        options += " -mask_data %s" % (mask_data if isinstance(mask_data, str) else ",".join(mask_data)) \
            if mask_data is not None else ""

        self.execute(options)

    def make_protein_db(self, input_file, db_title, mask_data, output_file=None, input_format="fasta", parse_seqids=True):
        # mask_data can be either string or list of strings
        self.make_db(input_file, db_title, mask_data, "prot", output_file=output_file,
                     input_format=input_format, parse_seqids=parse_seqids)

    def make_nucleotide_db(self, input_file, db_title, mask_data, output_file=None, input_format="fasta", parse_seqids=True):
        # mask_data can be either string or list of strings
        self.make_db(input_file, db_title, mask_data, "nucl", output_file=output_file,
                     input_format=input_format, parse_seqids=parse_seqids)


class BLASTDbCmd(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blastdbcmd", path=path, max_threads=max_threads)

    def db_info(self, db_name):
        options = " -info"
        options += " -db %s" % db_name

        self.execute(options)


if __name__ == "__main__":

    blast_plus = BLASTPlus()
    workdir = "/home/mahajrod/genetics/nxf/annotation/test/Dmel/nxf1/"
    os.chdir(workdir)
    blast_plus.make_blast_plus_db("/home/mahajrod/genetics/nxf/annotation/test/Dmel/nxf1/dmel-all-chromosome-r5.54.fasta",
                                  "dmel_all_chromosome.asnb", "dmel_all_chromosome")