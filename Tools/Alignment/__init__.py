__author__ = 'mahajrod'

from Tools.Alignment import Bowtie2, BWA

max_threads = 4
bowtie2_path = ""
bwa_path = ""

Bowtie2 = Bowtie2(path=bowtie2_path, max_threads=max_threads)
BWA = BWA(path=bwa_path, max_threads=max_threads)