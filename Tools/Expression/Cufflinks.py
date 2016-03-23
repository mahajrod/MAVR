#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class Gffread(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "gffread", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(input_gff_file, genomic_fasta_file=None, output_cds_file=None,
                             output_transcripts_file=None):

        options = " %s" % input_gff_file
        options += " -g %s" % genomic_fasta_file if genomic_fasta_file else ""
        options += " -x %s" % output_cds_file if output_cds_file else ""
        options += " -w %s" % output_transcripts_file if output_transcripts_file else ""

        return options

    def extract_transcript_sequences(self, input_gff_file, genomic_fasta_file, output_prefix):

        output_cds_file = "%s.cds" % output_prefix
        output_transcripts_file = "%s.transcript" % output_prefix
        options = self.parse_common_options(input_gff_file,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_cds_file=output_cds_file,
                                            output_transcripts_file=output_transcripts_file)
        #print options
        self.execute(options=options)

if __name__ == "__main__":
    pass
