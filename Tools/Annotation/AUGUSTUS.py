#!/usr/bin/env python
import os
import shutil

from Routines import FileRoutines, SequenceRoutines, MatplotlibRoutines
from CustomCollections.GeneralCollections import IdList, SynDict

from Tools.Abstract import Tool
from Tools.LinuxTools import CGAS


class AUGUSTUS(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "augustus", path=path, max_threads=max_threads)

    def parse_options(self, species, genome_file="", strand="both", gene_model=None, output_gff3=True,
                      other_options="", config_dir=None, use_softmasking=None, hints_file=None,
                      extrinsicCfgFile=None, predict_UTR=None):

        """
            AUGUSTUS (3.0.2) is a gene prediction tool for eukaryotes
        written by Mario Stanke (mstanke@gwdg.de) and Oliver Keller (keller@cs.uni-goettingen.de).

        usage:
        augustus [parameters] --species=SPECIES queryfilename

        'queryfilename' is the filename (including relative path) to the file containing the query sequence(s)
        in fasta format.

        SPECIES is an identifier for the species. Use --species=help to see a list.

        parameters:
        --strand=both, --strand=forward or --strand=backward
        --genemodel=partial, --genemodel=intronless, --genemodel=complete, --genemodel=atleastone or --genemodel=exactlyone
          partial      : allow prediction of incomplete genes at the sequence boundaries (default)
          intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes
          complete     : only predict complete genes
          atleastone   : predict at least one complete gene
          exactlyone   : predict exactly one complete gene
        --singlestrand=true
          predict genes independently on each strand, allow overlapping genes on opposite strands
          This option is turned off by default.
        --hintsfile=hintsfilename
          When this option is used the prediction considering hints (extrinsic information) is turned on.
          hintsfilename contains the hints in gff format.
        --AUGUSTUS_CONFIG_PATH=path
          path to config directory (if not specified as environment variable)
        --alternatives-from-evidence=true/false
          report alternative transcripts when they are suggested by hints
        --alternatives-from-sampling=true/false
          report alternative transcripts generated through probabilistic sampling
        --sample=n
        --minexonintronprob=p
        --minmeanexonintronprob=p
        --maxtracks=n
          For a description of these parameters see section 4 of README.TXT.
        --proteinprofile=filename
          When this option is used the prediction will consider the protein profile provided as parameter.
          The protein profile extension is described in section 7 of README.TXT.
        --progress=true
          show a progressmeter
        --gff3=on/off
          output in gff3 format
        --predictionStart=A, --predictionEnd=B
          A and B define the range of the sequence for which predictions should be found.
        --UTR=on/off
          predict the untranslated regions in addition to the coding sequence. This currently works only for a subset of species.
        --noInFrameStop=true/false
          Do not report transcripts with in-frame stop codons. Otherwise, intron-spanning stop codons could occur. Default: false
        --noprediction=true/false
          If true and input is in genbank format, no prediction is made. Useful for getting the annotated protein sequences.
        --uniqueGeneId=true/false
          If true, output gene identifyers like this: seqname.gN

        For a complete list of parameters, type "augustus --paramlist".
        An exhaustive description can be found in the file README.TXT.
            """

        options = " --uniqueGeneId=true"
        options += " --extrinsicCfgFile=%s" % extrinsicCfgFile if extrinsicCfgFile else ""
        options += " --softmasking=1" if use_softmasking else ""
        options += " --UTR=on" if predict_UTR else ""
        options += " --hintsfile=%s" % hints_file if hints_file else ""
        options += (" %s" % other_options) if other_options else ""
        options += " --strand=%s" % strand
        options += (" --genemodel=%s" % gene_model) if gene_model else ""
        options += " --gff3=%s" % ("on" if output_gff3 else "no")
        options += " --species=%s" % species
        options += (" --AUGUSTUS_CONFIG_PATH=%s" % config_dir) if config_dir else ""
        options += " %s" % genome_file

        return options

    def predict(self, species, genome_file, output, strand="both", gene_model="complete", output_gff3=True,
                other_options="", use_softmasking=None, hints_file=None, extrinsicCfgFile=None,
                predict_UTR=None):
        options = self.parse_options(species, genome_file=genome_file, strand=strand, gene_model=gene_model,
                                     output_gff3=output_gff3, other_options=other_options,
                                     use_softmasking=use_softmasking, hints_file=hints_file,
                                     extrinsicCfgFile=extrinsicCfgFile, predict_UTR=predict_UTR)
        options += " > %s" % output
        self.execute(options)

    def parallel_predict(self, species, genome_file, output, strand="both", gene_model=None, output_gff3=True,
                         other_options="", split_dir="splited_input", splited_output_dir="splited_output_dir",
                         config_dir=None, combine_output_to_single_file=True, use_softmasking=None, hints_file=None,
                         extrinsicCfgFile=None, predict_UTR=None):
        common_options = self.parse_options(species, genome_file="", strand=strand, gene_model=gene_model,
                                            output_gff3=output_gff3, other_options=other_options,
                                            config_dir=config_dir, use_softmasking=use_softmasking,
                                            hints_file=hints_file, extrinsicCfgFile=extrinsicCfgFile,
                                            predict_UTR=predict_UTR)

        splited_dir = FileRoutines.check_path(split_dir)
        splited_out_dir = FileRoutines.check_path(splited_output_dir)
        FileRoutines.save_mkdir(splited_dir)
        FileRoutines.save_mkdir(splited_out_dir)

        self.split_fasta_by_seq_len(genome_file, splited_dir)

        input_list_of_files = sorted(os.listdir(splited_dir))
        list_of_output_files = []
        options_list = []
        for filename in input_list_of_files:
            input_file = "%s%s" % (splited_dir, filename)
            output_file = "%s%s.gff" % (splited_out_dir, filename)
            list_of_output_files.append(output_file)
            options = common_options

            options += " %s" % input_file
            options += " > %s" % output_file
            options_list.append(options)

        self.parallel_execute(options_list)

        if combine_output_to_single_file:
            CGAS.cat(list_of_output_files, output=output)

    def extract_proteins_from_output(self, augustus_output, protein_output, evidence_stats_file=None,
                                     supported_by_hints_file=None, complete_proteins_id_file=None, id_prefix="p."):
        if evidence_stats_file:
            ev_fd = open(evidence_stats_file, "w")
            ev_fd.write("#gene_id\ttranscript_id\tsupported_fraction\tcds_support\tintron_support\t")
            ev_fd.write("5'UTR_support\t3'UTR_support\tincompatible_hints_groups\tprotein_length\n")

        if evidence_stats_file:
            sup_fd = open(supported_by_hints_file, "w")
            sup_fd.write("#gene_id\ttranscript_id\tsupported_fraction\tcds_support\tintron_support\t")
            sup_fd.write("5'UTR_support\t3'UTR_support\tincompatible_hints_groups\tprotein_length\n")

        if complete_proteins_id_file:
            complete_fd = open(complete_proteins_id_file, "w")

        with open(protein_output, "w") as out_fd:
            with open(augustus_output, "r") as in_fd:
                for line in in_fd:
                    if line[:12] == "# start gene":
                        gene = line.strip().split()[-1]
                    elif "\ttranscript\t" in line:
                        transcript_id = line.split("\t")[8].split(";")[0].split("=")[1]
                        start_presence = False
                        stop_presence = False
                        #out_fd.write(">%s%s\t gene=%s\n" % (id_prefix, transcript_id, gene))
                    elif "\tstart_codon\t" in line:
                        start_presence = True
                    elif "\tstart_codon\t" in line:
                        stop_presence = True
                    elif "# protein sequence" in line:
                        protein = line.strip().split("[")[-1]
                        if "]" in protein:
                            protein = protein.split("]")[0]
                        else:
                            while True:
                                part = in_fd.next().split()[-1]
                                if "]" in part:
                                    protein += part.split("]")[0]
                                    break
                                else:
                                    protein += part
                        if complete_proteins_id_file:
                            print "AAAAA"
                            if start_presence and stop_presence:
                                complete_fd.write("%s%s\n" % (id_prefix, transcript_id))

                        out_fd.write(">%s%s\t gene=%s\n start_presence=%s stop_presence=%s" % (id_prefix, transcript_id,
                                                                                               gene,
                                                                                               str(start_presence),
                                                                                               str(stop_presence)))
                        out_fd.write(protein)
                        protein_len = len(protein)
                        out_fd.write("\n")

                    elif evidence_stats_file or supported_by_hints_file:
                        if line[:17] == "# % of transcript":
                            supported_fraction = line.strip().split()[-1]
                            while True:
                                tmp_line = in_fd.next()
                                if tmp_line[:12] == "# CDS exons:":
                                    cds_support = tmp_line.strip().split()[-1]
                                elif tmp_line[:14] == "# CDS introns:":
                                    introns_support = tmp_line.strip().split()[-1]
                                elif tmp_line[:13] == "# 5'UTR exons":
                                    five_utr_support = tmp_line.strip().split()[-1]
                                elif tmp_line[:13] == "# 3'UTR exons":
                                    three_introns_support = tmp_line.strip().split()[-1]
                                elif tmp_line[:27] == "# incompatible hint groups:":
                                    incompatible_hint_groups = tmp_line.strip().split()[-1]
                                    if evidence_stats_file:
                                        ev_fd.write("%s\t%s\t%s\t" % (gene, transcript_id, supported_fraction))
                                        ev_fd.write("%s\t%s\t%s\t%s\t%s\t%i\n" % (cds_support, introns_support,
                                                                                  five_utr_support,
                                                                                  three_introns_support,
                                                                                  incompatible_hint_groups,
                                                                                  protein_len))
                                    if supported_by_hints_file and (float(supported_fraction) > 0):
                                        sup_fd.write("%s\t%s\t%s\t" % (gene, transcript_id, supported_fraction))
                                        sup_fd.write("%s\t%s\t%s\t%s\t%s\t%i\n" % (cds_support, introns_support,
                                                                                   five_utr_support,
                                                                                   three_introns_support,
                                                                                   incompatible_hint_groups,
                                                                                   protein_len))

                                    break

        if evidence_stats_file:
            ev_fd.close()

        self.extract_longest_isoforms(evidence_stats_file, "%s.longest_pep" % evidence_stats_file,
                                      minimum_supported_fraction=0)
        SequenceRoutines.extract_sequence_by_ids(protein_output,
                                                 "%s.longest_pep.ids" % evidence_stats_file,
                                                 "%s.longest_pep.pep" % evidence_stats_file)

        if supported_by_hints_file:
            supported_by_hints_longest_pep_evidence = "%s.longest_pep" % supported_by_hints_file
            supported_by_hints_longest_pep = "%s.longest_pep.pep" % supported_by_hints_file
            supported_by_hints_longest_pep_ids = "%s.longest_pep.ids" % supported_by_hints_file
            self.extract_longest_isoforms(evidence_stats_file, supported_by_hints_longest_pep_evidence,
                                          minimum_supported_fraction=0.00001)
            SequenceRoutines.extract_sequence_by_ids(protein_output,
                                                     supported_by_hints_longest_pep_ids,
                                                     supported_by_hints_longest_pep)

        evidence_files = (evidence_stats_file,
                          "%s.longest_pep" % evidence_stats_file,
                          "%s.longest_pep" % supported_by_hints_file) if supported_by_hints_file else \
                          (evidence_stats_file,)
        for evidence_file in evidence_files:
            print ("Drawing transcript support distribution for %s" % evidence_file)
            MatplotlibRoutines.percent_histogram_from_file(evidence_file,
                                                           evidence_file,
                                                           column_list=(2,), separator=None,
                                                           comments="#", n_bins=20,
                                                           title="Transcript support by hints",
                                                           xlabel="%%", ylabel="Number",
                                                           extensions=["svg", "png"],
                                                           legend_location="upper center",
                                                           stats_as_legend=True)

    @staticmethod
    def extract_evidence_by_ids(evidence_file, id_file, output_evidence_file, mode="transcript"):
        # possible modes: transcript, gene
        ids = IdList()
        ids.read(id_file, comments_prefix="#")

        column_id = 0 if mode == "gene" else 1

        with open(evidence_file, "r") as ev_fd:
            with open(output_evidence_file, "w") as out_fd:
                for line in ev_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                        continue

                    entry_id = line.split("\t")[column_id]
                    if entry_id in ids:
                        out_fd.write(line)

    @staticmethod
    def extract_longest_isoforms(evidence_file, filtered_evidence_file, minimum_supported_fraction=0):
        longest_id_file = "%s.ids" % filtered_evidence_file
        ev_fd = open(evidence_file, "r")
        with open(longest_id_file, "w") as id_fd:
                with open(filtered_evidence_file, "w") as filtered_ev_fd:
                    filtered_ev_fd.write(ev_fd.readline())
                    for line in ev_fd:
                        prev_line = line
                        line_list = prev_line.strip().split("\t")
                        if float(line_list[2]) < minimum_supported_fraction:
                            continue
                        prev_gene = line_list[0]
                        prev_transcript = line_list[1]
                        prev_pep_len = line_list[-1]

                        break
                    for line in ev_fd:

                        line_list = line.strip().split("\t")
                        if float(line_list[2]) < minimum_supported_fraction:
                            continue
                        gene = line_list[0]
                        transcript = line_list[1]
                        pep_len = line_list[-1]
                        if gene == prev_gene:
                            if pep_len > prev_pep_len:
                                prev_line = line
                                prev_gene = gene
                                prev_transcript = transcript
                                prev_pep_len = pep_len
                            continue
                        id_fd.write("%s\n" % prev_transcript)
                        filtered_ev_fd.write(prev_line)
                        prev_line = line
                        prev_gene = gene
                        prev_transcript = transcript
                        prev_pep_len = pep_len
                    # if script breaks here it means that prediction was made without hints or all genes doesnt have hint support
                    id_fd.write("%s\n" % prev_transcript)
        ev_fd.close()

    @staticmethod
    def extract_CDS_annotations_from_output(augustus_output, CDS_output):
        CGAS.grep("'\\tCDS\\t'", augustus_output, output=CDS_output, use_regexp=True)

    @staticmethod
    def extract_gene_ids_from_output(augustus_output, all_genes_output):
        CGAS.cgas(augustus_output, grep_pattern="'\\tgene\\t'", sed_string="'s/.*ID=//'", output=all_genes_output,
                  grep_use_regexp=True)

    def exonerate_to_hints(self, exonerate_gff, hint_gff, priority=None, min_intron_len=None, max_intron_len=None,
                           CDS_part_cutoff=None, source=None):

        options = " --in=%s" % exonerate_gff
        options += " --out=%s" % hint_gff
        options += " --priority=%i" % priority if priority else ""
        options += " --minintronlen=%i" % min_intron_len if min_intron_len else ""
        options += " --maxintronlen=%i" % max_intron_len if max_intron_len else ""
        options += " --CDSpart_cutoff=%i" % CDS_part_cutoff if CDS_part_cutoff else ""
        options += " --source=%s" % source if source else ""

        self.execute(options, cmd="exonerate2hints.pl")

    def join_multiple_hints(self, in_file, out_file):
        options = " > %s" % out_file
        # sorting of hints is necessary before merging
        cmd = "cat %s | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl" % in_file
        self.execute(options, cmd=cmd)

    def assign_synonyms_to_features_from_augustus_gff(self, input_gff, output_file, prefix,
                                                      number_of_digits_in_number=8, feature_type="gene"):
        options = " > %s" % output_file
        cmd = """grep -P "\\t%s\\t" %s | sed 's/.*ID=//;s/;.*//' | awk -F'\\t' 'BEGIN {NUMBER=1};{printf "%%s\\t%s%%0%ii\\n",$1,NUMBER; NUMBER=NUMBER+1;}'""" % (feature_type, input_gff, prefix, number_of_digits_in_number)

        self.execute(options, cmd=cmd)

    def assign_synonyms_to_annotations_from_augustus_gff(self, input_gff, output_prefix, species_prefix,
                                                         number_of_digits_in_number=8):

        for feature in "gene", "transcript":
            output_file = "%s.%s.syn" % (output_prefix, feature)
            feature_prefix = "%s%s" % (species_prefix.upper(), "G" if feature == "gene" else "T")
            options = " > %s" % output_file
            cmd = """grep -P "\\t%s\\t" %s | sed 's/.*ID=//;s/;.*//' | awk -F'\\t' 'BEGIN {NUMBER=1};{printf "%%s\\t%s%%0%ii\\n",$1,NUMBER; NUMBER=NUMBER+1;}'""" % (feature, input_gff, feature_prefix, number_of_digits_in_number)
            self.execute(options, cmd=cmd)

        options = " %s.transcript.syn > %s.protein.syn" % (output_prefix, output_prefix)

        cmd = "sed 's/%sT/%sP/'" % (species_prefix, species_prefix)

        self.execute(options, cmd=cmd)

        options = " %s.transcript.syn > %s.cds.syn" % (output_prefix, output_prefix)
        cmd = "sed 's/%sT/%sC/;s/\\t/.cds\\t/'" % (species_prefix, species_prefix)
        self.execute(options, cmd=cmd)

    @staticmethod
    def replace_augustus_ids_by_syn(augustus_gff, output_gff, genes_syn_file, transcripts_syn_file,
                                    cds_syn_file=None):

        genes_syn_dict = SynDict()
        genes_syn_dict.read(genes_syn_file, comments_prefix="#")
        transcripts_syn_dict = SynDict()
        transcripts_syn_dict.read(transcripts_syn_file, comments_prefix="#")
        cds_syn_dict = SynDict()
        if cds_syn_file:
            cds_syn_dict.read(cds_syn_file, comments_prefix="#")
        with open(augustus_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    tmp = line.strip()
                    if len(tmp) < 13:
                        out_fd.write(line)
                        continue
                    if tmp[:12] != "# start gene":
                        out_fd.write(line)
                        continue
                    augustus_gene_id = tmp.split(" ")[-1]
                    gene_syn_id = genes_syn_dict[augustus_gene_id]
                    augustus_transcript_id = ""
                    augustus_transcript_parent = ""
                    out_fd.write("# start gene %s\n" % gene_syn_id)
                    tmp = in_fd.next().strip()
                    while True:
                        while tmp[0] != "#":
                            tmp_list = tmp.split("\t")
                            feature_type = tmp_list[2]
                            edited_str = "\t".join(tmp_list[:-1])
                            info_field_list = tmp_list[-1].split(";")
                            if feature_type == "gene":
                                edited_str += "\tID=%s\n" % gene_syn_id
                            elif feature_type == "transcript":
                                for entry in info_field_list:
                                    if "ID" in entry:
                                        augustus_transcript_id = entry.split("=")[-1]
                                        transcript_syn_id = transcripts_syn_dict[augustus_transcript_id]
                                    if "Parent" in entry:
                                        augustus_transcript_parent = entry.split("=")[-1]
                                        if augustus_transcript_parent != augustus_gene_id:
                                            raise ValueError("Transcript parent id and gene id are not same!")
                                edited_str += "\tID=%s;Parent=%s\n" % (transcript_syn_id, gene_syn_id)
                            elif feature_type == "CDS":
                                for entry in info_field_list:
                                    if "ID" in entry:
                                        augustus_cds_id = entry.split("=")[-1]
                                        cds_syn_id = cds_syn_dict[augustus_cds_id] if cds_syn_dict else "%s.cds" % transcripts_syn_dict[augustus_cds_id[:-4]]
                                    if "Parent" in entry:
                                        augustus_cds_parent = entry.split("=")[-1]
                                        if augustus_cds_parent != augustus_transcript_id:
                                            raise ValueError("CDS parent id and transcript id are not same!")
                                edited_str += "\tID=%s;Parent=%s\n" % (cds_syn_id, transcript_syn_id)
                            elif (feature_type == "stop_codon") or (feature_type == "start_codon"):
                                for entry in info_field_list:
                                    if "Parent" in entry:
                                        augustus_feature_parent = entry.split("=")[-1]
                                        if augustus_feature_parent != augustus_transcript_id:
                                            raise ValueError("Feature parent id and transcript id are not same!")
                                edited_str += "\tParent=%s\n" % (transcript_syn_id)
                            else:
                                edited_str = tmp

                            out_fd.write(edited_str)
                            tmp = in_fd.next().strip()
                        while tmp[0] == "#":
                            if "# end gene" in tmp:
                                break
                            out_fd.write(tmp + "\n")
                            tmp = in_fd.next().strip()
                        if "# end gene" in tmp:
                                break
                    out_fd.write("# end gene %s\n" % gene_syn_id)

    @staticmethod
    def replace_augustus_ids(augustus_gff, output_prefix, species_prefix=None, number_of_digits_in_id=8):

        output_gff = "%s.renamed.gff" % output_prefix
        genes_syn_file = "%s.gene.syn" % output_prefix
        transcripts_syn_file = "%s.transcript.syn" % output_prefix
        cds_syn_file = "%s.cds.syn" % output_prefix
        genes_syn_dict = SynDict()
        transcripts_syn_dict = SynDict()
        cds_syn_dict = SynDict()
        gene_counter = 0
        gene_id_template = "%sG%%0%ii" % (species_prefix, number_of_digits_in_id)
        transcripts_counter = 0
        transcript_id_template = "%sT%%0%ii" % (species_prefix, number_of_digits_in_id)
        cds_counter = 0
        cds_id_template = "%sC%%0%ii" % (species_prefix, number_of_digits_in_id)
        with open(augustus_gff, "r") as in_fd:
            with open(output_gff, "w") as out_fd:
                for line in in_fd:
                    tmp = line.strip()
                    if len(tmp) < 13:
                        out_fd.write(line)
                        continue
                    if tmp[:12] != "# start gene":
                        out_fd.write(line)
                        continue
                    augustus_gene_id = tmp.split(" ")[-1]
                    gene_counter += 1

                    gene_syn_id = gene_id_template % gene_counter
                    genes_syn_dict[augustus_gene_id] = gene_syn_id
                    augustus_transcript_id = ""
                    augustus_transcript_parent = ""
                    out_fd.write("# start gene %s\n" % gene_syn_id)
                    tmp = in_fd.next().strip()
                    while True:
                        while tmp[0] != "#":
                            tmp_list = tmp.split("\t")
                            feature_type = tmp_list[2]
                            edited_str = "\t".join(tmp_list[:-1])
                            info_field_list = tmp_list[-1].split(";")
                            if feature_type == "gene":
                                edited_str += "\tID=%s\n" % gene_syn_id
                            elif feature_type == "transcript":
                                for entry in info_field_list:
                                    if "ID" in entry:
                                        augustus_transcript_id = entry.split("=")[-1]
                                        if augustus_transcript_id not in transcripts_syn_dict:
                                            transcripts_counter += 1
                                            transcripts_syn_dict[augustus_transcript_id] = transcript_id_template % transcripts_counter
                                        transcript_syn_id = transcripts_syn_dict[augustus_transcript_id]
                                    if "Parent" in entry:
                                        augustus_transcript_parent = entry.split("=")[-1]
                                        if augustus_transcript_parent != augustus_gene_id:
                                            raise ValueError("Transcript parent id and gene id are not same!")
                                edited_str += "\tID=%s;Parent=%s\n" % (transcript_syn_id, gene_syn_id)
                            elif feature_type == "CDS":
                                for entry in info_field_list:
                                    if "ID" in entry:
                                        augustus_cds_id = entry.split("=")[-1]
                                        if augustus_cds_id not in cds_syn_dict:
                                            cds_counter += 1
                                            cds_syn_dict[augustus_cds_id] = cds_id_template % cds_counter
                                        cds_syn_id = cds_syn_dict[augustus_cds_id]
                                    if "Parent" in entry:
                                        augustus_cds_parent = entry.split("=")[-1]
                                        if augustus_cds_parent != augustus_transcript_id:
                                            raise ValueError("CDS parent id and transcript id are not same!")
                                edited_str += "\tID=%s;Parent=%s\n" % (cds_syn_id, transcript_syn_id)
                            elif (feature_type == "stop_codon") or (feature_type == "start_codon"):
                                for entry in info_field_list:
                                    if "Parent" in entry:
                                        augustus_feature_parent = entry.split("=")[-1]
                                        if augustus_feature_parent != augustus_transcript_id:
                                            raise ValueError("Feature parent id and transcript id are not same!")
                                edited_str += "\tParent=%s\n" % transcript_syn_id
                            else:
                                edited_str = tmp

                            out_fd.write(edited_str)
                            tmp = in_fd.next().strip()
                        while tmp[0] == "#":
                            if "# end gene" in tmp:
                                break
                            out_fd.write(tmp + "\n")
                            tmp = in_fd.next().strip()
                        if "# end gene" in tmp:
                                break
                    out_fd.write("# end gene %s\n" % gene_syn_id)
        genes_syn_dict.write(genes_syn_file)
        transcripts_syn_dict.write(transcripts_syn_file)
        cds_syn_dict.write(cds_syn_file)





