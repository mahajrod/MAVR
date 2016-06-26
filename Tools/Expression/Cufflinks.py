#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class Gffread(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "gffread", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(input_gff_file, genomic_fasta_file=None, output_cds_file=None, output_protein_file=None,
                             output_transcripts_file=None, coding_only=False):

        options = " %s" % input_gff_file
        options += " -C" if coding_only else ""
        options += " -g %s" % genomic_fasta_file if genomic_fasta_file else ""
        options += " -x %s" % output_cds_file if output_cds_file else ""
        options += " -w %s" % output_transcripts_file if output_transcripts_file else ""
        options += " -y %s" % output_protein_file if output_protein_file else ""
        return options

    def extract_transcript_sequences(self, input_gff_file, genomic_fasta_file, output_prefix, coding_only=False):

        output_cds_file = "%s.cds" % output_prefix
        output_transcripts_file = "%s.transcript" % output_prefix
        output_protein_file = "%s.protein" % output_prefix
        options = self.parse_common_options(input_gff_file,
                                            coding_only=coding_only,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_cds_file=output_cds_file,
                                            output_protein_file=output_protein_file,
                                            output_transcripts_file=output_transcripts_file)
        #print options
        self.execute(options=options)

    def extract_cds(self, input_gff_file, genomic_fasta_file, output_file):

        options = self.parse_common_options(input_gff_file,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_cds_file=output_file,
                                            coding_only=True,
                                            output_transcripts_file=None)
        #print options
        self.execute(options=options)

    def extract_mRNA(self, input_gff_file, genomic_fasta_file, output_file, coding_only=False):

        options = self.parse_common_options(input_gff_file,
                                            coding_only=coding_only,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_transcripts_file=output_file)
        #print options
        self.execute(options=options)

if __name__ == "__main__":
    pass
