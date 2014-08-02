#!/usr/bin/env python
import os


class Bowtie2():
    def __init__(self):
        pass

    def index(self, reference, index_name):
        os.system("bowtie2-build %s %s" % (reference, index_name))

    def align(self,
              bowtie2_index,
              forward_reads,
              reverse_reads=None,
              max_threads=5,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              output_file="alignment.sam"):

        if reverse_reads:
            reads = "-1 %s -2 %s" % (forward_reads, reverse_reads)
        else:
            reads = "-U %s" % forward_reads

        os.system("bowtie2 --%s --%s -p %i -x %s %s > %s"
                  % (alignment_mode, quality_score, max_threads, bowtie2_index, reads, output_file))


class BWA():
    def __init__(self):
        pass

    def make_index(self, fasta_file, prefix=None):
        index_prefix = ""
        if prefix:
            index_prefix = "-p %s" % prefix
        os.system("bwa index %s %s" % (index_prefix, fasta_file))

    def align_mem(self, index, forward_reads, reverse_reads=None, output_file="alignment.sam", max_threads=4):
        reads = forward_reads
        if reverse_reads:
            reads += " %s" % reverse_reads
        os.system("bwa mem -t %i %s %s > %s" % (max_threads, index, reads, output_file))


if __name__ == "__main__":
    os.chdir("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference")
    bwa = BWA()
    bwa.make_index("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24s.fasta",
                   prefix="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24s")