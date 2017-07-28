__author__ = 'mahajrod'

from Bio import AlignIO

from Bio.Phylo.Consensus import bootstrap

from Routines.MultipleAlignment import MultipleAlignmentRoutines


class PhylogeneticsRoutines(MultipleAlignmentRoutines):
    def __init__(self):
        MultipleAlignmentRoutines.__init__(self)

    @staticmethod
    def bootstrap_alignment(alignment_file, output_directory, output_prefix,
                            replicate_number, format="fasta"):
        alignment = AlignIO.read(alignment_file, format=format)

        bootstrap_replicates = bootstrap(alignment, replicate_number)

        for i in range(0, len(bootstrap_replicates)):
            output_file = "%s/%s_%i.%s" % (output_directory, output_prefix, i+1, format)
            AlignIO.write(bootstrap_replicates[i], output_file, format=format)
