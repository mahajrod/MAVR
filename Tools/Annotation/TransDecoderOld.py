#!/usr/bin/env python
import os

from Routines.File import split_filename

from Tools.Abstract import Tool


class TransDecoderOld(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "TransDecoder", path=path, max_threads=max_threads)

    def extract_longest_orfs(self, transcripts_file, genetic_code=None, analyze_only_top_strand=None,
                             minimum_protein_length=100):
        """
        Possible genetic codes:
            universal (default)
            Euplotes
            Tetrahymena
            Candida
            Acetabularia
            Mitochondrial-Canonical
            Mitochondrial-Vertebrates
            Mitochondrial-Arthropods
            Mitochondrial-Echinoderms
            Mitochondrial-Molluscs
            Mitochondrial-Ascidians
            Mitochondrial-Nematodes
            Mitochondrial-Platyhelminths
            Mitochondrial-Yeasts
            Mitochondrial-Euascomycetes
            Mitochondrial-Protozoans
        """
        options = " -t %s" % transcripts_file
        options += " -G %s" % genetic_code if genetic_code else ""
        options += " -S" if analyze_only_top_strand else ""
        options += " -m %i" % minimum_protein_length

        self.execute(options, cmd="TransDecoder.LongOrfs")

    def predict_pep(self, transcripts_file, minimum_orf_length_if_no_other_evidence=None,
                    pfam_hits=None, blastp_hits=None,
                    file_with_orfs_for_training=None,
                    number_of_top_orfs_for_training=None):

        options = " -t %s" % transcripts_file
        options += " --retain_long_orfs %i" % minimum_orf_length_if_no_other_evidence
        options += " --retain_pfam_hits %s" % pfam_hits if pfam_hits else None
        options += " --retain_blastp_hits %s" % blastp_hits if blastp_hits else None
        options += " --train %s" % file_with_orfs_for_training if file_with_orfs_for_training else None
        options += " -T %i" % number_of_top_orfs_for_training if number_of_top_orfs_for_training else None

        self.execute(options, cmd="TransDecoder.Predict")




    """
    def annotate_pep(self, transcripts_file, genetic_code=None, analyze_only_top_strand=None,
                     reuse_existing_files=None, minimum_protein_length=100, write_temp__to_current_dir=None):


        options = " --CPU %i" % self.threads
        options += " -t %s" % transcripts_file

        options += " -G %s" % genetic_code if genetic_code else ""
        options += " -S" if analyze_only_top_strand else ""
        options += " --reuse" if reuse_existing_files else ""
        options += " -m %i" % minimum_protein_length

        options += " --workdir" if write_temp__to_current_dir else ""
        options += " > %s" % output

        self.execute(options)

    """