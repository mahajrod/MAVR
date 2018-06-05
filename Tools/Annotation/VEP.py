#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool
from Parsers.Barrnap import CollectionBARRNAP


class VEP(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "vep", path=path, max_threads=max_threads)

    def parse_options(self):
        """
        if kingdom not in self.kingdoms:
            raise ValueError("Wrong kingdom of life. Allowed: euk, bac, mito, arc")

        options = " --threads %i" % self.threads
        options += " --kingdom %s" % kingdom
        options += " --lencutoff %f" % length_cutoff if length_cutoff else ""
        options += " --reject %f" % reject_cutoff if reject_cutoff else ""
        options += " --evalue %f" % evalue_cutoff if evalue_cutoff else ""
        options += " %s" % input_fasta
        options += " > %s " % gff_file
        options += " 2> %s" % log_file

        return options
        """
        pass

    def prepare_gff(self, input_gff_list, output_prefix):

        sorted_gff = "%s.sorted.gff" % output_prefix
        bgzip_compressed_gff = "%s.sorted.gff.bz" % output_prefix
        sorting_string = "cat %s | grep -v '#' | sort -k1,1 -k4,4n -k5,5n | tee %s | bgzip -c > %s" % (input_gff_list if isinstance(input_gff_list, str) else " ".join(input_gff_list),
                                                                                                                sorted_gff,
                                                                                                                bgzip_compressed_gff)
        indexing_string = "tabix -p gff %s" % bgzip_compressed_gff
        self.execute(cmd=sorting_string)
        self.execute(cmd=indexing_string)

if __name__ == "__main__":
    pass
