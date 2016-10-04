#!/usr/bin/env python
from Tools.Abstract import Tool


class MaskFasta(Tool):
    """
    Tool:    bedtools getfasta (aka fastaFromBed)
    Version: v2.17.0
    Summary: Extract DNA sequences into a fasta file based on feature coordinates.

    Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>

    Options:
        -fi	Input FASTA file
        -bed	BED/GFF/VCF file of ranges to extract from -fi
        -fo	Output file (can be FASTA or TAB-delimited)
        -name	Use the name field for the FASTA header
        -split	given BED12 fmt., extract and concatenate the sequencesfrom the BED "blocks" (e.g., exons)
        -tab	Write output in TAB delimited format.
            - Default is FASTA format.

        -s	Force strandedness. If the feature occupies the antisense,
            strand, the sequence will be reverse complemented.
            - By default, strand information is ignored.

    """

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools maskfasta", path=path, max_threads=max_threads)

    def mask(self, input_fasta, output_fasta, masking_gff, softmasking=False, masking_character=None):

        options = " -bed %s" % masking_gff
        options += " -fi %s" % input_fasta
        options += " -fo %s" % output_fasta
        options += " -soft" if softmasking else ""
        options += (" -mc %s" % masking_character) if masking_character else ""

        self.execute(options)
