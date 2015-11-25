#!/usr/bin/env python
import os

from Routines.File import split_filename, check_path, save_mkdir

from Tools.Abstract import Tool
from Tools.LinuxTools import CGAS


class AUGUSTUS(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "augustus", path=path, max_threads=max_threads)

    def parse_options(self, species, genome_file="", strand="both", gene_model=None, output_gff3=True,
                      other_options="", config_dir=None):

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
        options += (" %s" % other_options) if other_options else ""
        options += " --strand=%s" % strand
        options += (" --genemodel=%s" % gene_model) if gene_model else ""
        options += " --gff3=%s" % ("on" if output_gff3 else "no")
        options += " --species=%s" % species
        options += (" --AUGUSTUS_CONFIG_PATH=%s" % config_dir) if config_dir else ""
        options += " %s" % genome_file

        return options

    def predict(self, species, genome_file, output, strand="both", gene_model="complete", output_gff3=True,
                other_options=""):
        options = self.parse_options(species, genome_file=genome_file, strand=strand, gene_model=gene_model,
                                     output_gff3=output_gff3, other_options=other_options)
        options += " > %s" % output
        self.execute(options)

    def parallel_predict(self, species, genome_file, output, strand="both", gene_model=None, output_gff3=True,
                         other_options="", split_dir="splited_input", splited_output_dir="splited_output_dir",
                         config_dir=None, combine_output_to_single_file=True):
        common_options = self.parse_options(species, genome_file="", strand=strand, gene_model=gene_model,
                                            output_gff3=output_gff3, other_options=other_options,
                                            config_dir=config_dir)

        splited_dir = check_path(split_dir)
        splited_out_dir = check_path(splited_output_dir)
        save_mkdir(splited_dir)
        save_mkdir(splited_out_dir)

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

    @staticmethod
    def extract_proteins_from_output(augustus_output, protein_output, id_prefix="p."):
        with open(protein_output, "w") as out_fd:
            with open(augustus_output, "r") as in_fd:
                for line in in_fd:
                    if line[:12] == "# start gene":
                        gene = line.strip().split()[-1]
                        out_fd.write(">%s%s\t gene=%s\n" % (id_prefix, gene, gene))
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

                        out_fd.write(protein)
                        out_fd.write("\n")

    @staticmethod
    def extract_CDS_annotations_from_output(augustus_output, CDS_output):
        CGAS.grep("'\\tCDS\\t'", augustus_output, output=CDS_output, use_regexp=True)






