#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class LRNAScaffolder(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "SSPACE_Standard_v3.0.pl", path=path, max_threads=max_threads)

    def scaffold(self, genome_fasta, psl_file, scripts_dir=None, max_hit_length_percentage_for_single_hit=None):
        """
        IMPORTANT: Extension of contigs during scaffolding doesnt work in SSPACE 3.0!!!!!!!!!!!!!

        """

        options = " -d %s" % scripts_dir if scripts_dir else self.path
        options += " -i %s" % psl_file
        options += " -j %s" % genome_fasta
        options += " -l %i" % if extend_contigs else ""
        options += " -m %i" % min_overlap_len if min_overlap_len else ""
        options += " -o %i" % min_reads_for_contig_extension if min_reads_for_contig_extension else ""
        options += " -z %i" % min_contig_len if min_contig_len else ""
        options += " -k %i" % min_link_number if min_link_number else ""
        options += " -a %f" % min_link_ratio if min_link_ratio else ""
        options += " -n %i" % min_contig_overlap if min_contig_overlap else ""
        options += " -g %i" % max_gaps_for_bowtie if max_gaps_for_bowtie else ""
        options += " -S" if skip_read_processing else ""
        options += " -b %s" % output_basename if output_basename else ""
        options += " -v 1" if verbose_scaffolding else ""
        options += " -p 1" if make_dot_file else ""


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



        self.execute(options)

if __name__ == "__main__":
    pass
