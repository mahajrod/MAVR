__author__ = 'mahajrod'
from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC, Gapped


class MultipleAlignmentConverters:

    @staticmethod
    def convert_alignment(input_file, input_filetype, output_file, output_filetype, alphabet=Gapped(IUPAC.ambiguous_dna)):
        alignment = AlignIO.parse(input_file, input_filetype, alphabet=alphabet)
        AlignIO.write\
            (alignment, output_file, output_filetype)

    @staticmethod
    def fasta2paml(input, output):
        sequence_dict = SeqIO.to_dict(SeqIO.parse(input, format="fasta"))
        number_of_sequences = len(sequence_dict)
        alignment_length = len(sequence_dict[sequence_dict.keys()[0]])

        with open(output, "w") as out_fd:
            out_fd.write("%i %i\n" % (number_of_sequences, alignment_length))
            for record_id in sequence_dict:
                out_fd.write("%s\n%s\n" % (record_id, str(sequence_dict[record_id].seq)))