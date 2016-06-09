#!/usr/bin/env python
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

class MultipleAlignmentConverters:

    @staticmethod
    def convert_alignment(input_file, input_filetype, output_file, output_filetype, alphabet=Gapped(IUPAC.ambiguous_dna)):
        alignment = AlignIO.parse(input_file, input_filetype, alphabet=alphabet)
        AlignIO.write(alignment, output_file, output_filetype)