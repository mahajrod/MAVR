#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class SSPACE(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "SSPACE_Standard_v3.0.pl", path=path, max_threads=max_threads)

    def scaffold(self, library_file, contig_fasta, output_basename, extend_contigs=True, min_overlap_len=None,
                 min_reads_for_contig_extension=None, min_contig_len=None, min_link_number=None,
                 min_link_ratio=None, min_contig_overlap=None, max_gaps_for_bowtie=None, skip_read_processing=None,
                 verbose_scaffolding=None,
                 make_dot_file=True):

        options = " -T %i" % self.threads
        options += " -l %s" % library_file
        options += " -s %s" % contig_fasta
        options += " -x 1" if extend_contigs else ""
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

        self.execute(options)

if __name__ == "__main__":
    pass
