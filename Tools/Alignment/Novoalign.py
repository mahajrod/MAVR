#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class Novoalign(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "novoalign", path=path, max_threads=max_threads)

    def index(self, index_name, reference_files, kmer_size=None, step_size=None, soft_masking=False,
              bisulfite_index=False, color_index=False, internal_name=None):

        options = ""
        options += " -k %i" % kmer_size if kmer_size else ""
        options += " -s %i" % step_size if step_size else ""
        options += " -m" if soft_masking else ""
        options += " -b" if bisulfite_index else ""
        options += " -c" if color_index else ""
        options += " -n %s" % internal_name if internal_name else ""
        options += " %s" % index_name
        options += " %s" % " ".join(reference_files) if isinstance(reference_files, list) else " %s" % reference_files

        self.execute(options, cmd="novoindex")

    def align(self, reads, index, output_file="alignment.sam", gap_extend_pen=None,
              gap_open_pen=None, out_fmt="SAM", in_fmt="ILM1.8"):
        """
        -o %s       Output format : Native, Pairwise, SAM, extended
        -F %s       Input format:
                    FA              fasta without quality values
                    SLXFQ           fastq with solexa style quality values
                    STDFQ           fastq with Sanger style quality values
                    ILMFQ           fastq with Illumina-64 style quality values
                    ILM1.8          fastq with Illumina-33 style quality values
                    PRB
                    PRBnSEQ
                    QSEQ
                    BAMSE
                    BAMPE
                    BAM
        -g %i       Gap open penalty. Default 40.
        -x %i       Gap extend penalty. Default 6.
        -d %s       full path to index
        -f %s [%s]  Files with reads
        """
        options = ""
        options += " -g %i" % gap_open_pen if gap_open_pen else ""
        #options += " -x %i" % gap_extend_pen if gap_extend_pen else ""
        options += " -o %s" % out_fmt
        options += " -F %s" % in_fmt
        options += " -d %s" % index
        options += " -f %s" % " ".join(reference) if isinstance(reads, list) else " -f %s" % reads
        options += " -c %i" % self.max_threads
        options += " > %s" % output_file
        self.execute(options, cmd="novoalign")


if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/"
    reference = "Alsu24mc.fasta"
    index_name = "Alsu24mc.nidx"
    os.chdir(workdir)
    Novoalign = Novoalign()
    Novoalign.index(index_name, reference)
