__author__ = 'mahajrod'

from Bio import AlignIO

from Bio.Phylo.Consensus import bootstrap

from Routines.MultipleAlignment import MultipleAlignmentRoutines


class PhylogeneticsRoutines(MultipleAlignmentRoutines):
    def __init__(self):
        MultipleAlignmentRoutines.__init__(self)

    def bootstrap_alignment(self, alignment_file, output_directory, output_prefix,
                            replicate_number, format="fasta"):
        self.safe_mkdir(output_directory)
        alignment = AlignIO.read(alignment_file, format=format)

        bootstrap_replicate_generator = bootstrap(alignment, replicate_number)

        i = 0
        for replicate in bootstrap_replicate_generator:
            i += 1
            output_file = "%s/%s_%i.%s" % (output_directory, output_prefix, i, format)
            AlignIO.write(replicate, output_file, format=format)
