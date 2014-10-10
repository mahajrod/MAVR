#!/usr/bin/env python
from Tools.Abstract import Tool


class GetFasta(Tool):
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
        Tool.__init__(self, "bedtools getfasta", path=path, max_threads=max_threads)

    def get(self, bed_file, fasta_file, output_file, use_strand=False, use_name_field_for_fasta_header=False,
            output_tab=False):

        options = ""

        options += " -bed %s" % bed_file
        options += " -fi %s" % fasta_file
        options += " -fo %s" % output_file
        options += " -s" if use_strand else ""
        options += " -name" if use_name_field_for_fasta_header else ""
        options += " -tab" if output_tab else ""

        self.execute(options)

