#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class LRNAScaffolder(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "L_RNA_scaffolder.sh", path=path, max_threads=max_threads)

    def scaffold(self, genome_fasta, psl_file, output_dir="./", scripts_dir=None, max_hit_length_percentage_for_single_hit=None,
                 alignment_threshold_identity=None, max_intron_length=None, min_link_number=None):
        """
           -d        C++ and Perl programs directory path, -d is mandatory.
           -i        PSL file (result of BLAT), -i is mandatory.
           -j        Pre-assembly contig FASTA file (Database of BLAT), -j is mandatory.
           -l        length_coverage, the coverage threshold of alignment length to full length (Default is 0.95).
           -p        identity, the identity threshold of alignment (Default is 0.9).
           -o        output directory, where you will put assembly result, default is current directory,
                     you can mkdir an output-directory before assembly.
           -e        the maximal intron, default is 100000 bp.
           -f        frequency of the lowest routes, default is 1.
           -r        overlap contigs according to AGP file.
           -n        N number you want to indicate a gap if the gap is larger than the median intron size.


        """

        options = " -d %s" % scripts_dir if scripts_dir else self.path
        options += " -i %s" % psl_file
        options += " -j %s" % genome_fasta
        options += " -l %f" % max_hit_length_percentage_for_single_hit if max_hit_length_percentage_for_single_hit else "" #If one read has a hit of which length coverage was over the threshold, then this read would be filtered out.
        options += " -p %f" % alignment_threshold_identity if alignment_threshold_identity else "" #If one alignment has an identity over the threshold, then the alignment is kept for the further analysis.
        options += " -o %s" % output_dir if output_dir else ""
        options += " -e %i" % max_intron_length if max_intron_length else "" # default is 100000 bp.
        options += " -f %i" % min_link_number if min_link_number else ""    #default is 1.

        self.execute(options)

if __name__ == "__main__":
    pass
