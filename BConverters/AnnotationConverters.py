__author__ = 'mahajrod'
import os

from BCBio import GFF


class AnnotationConverters:

    @staticmethod
    def gff22gff3(input_file, output_file, target_lines=100000):

        in_fd = open(input_file, "r")
        out_fd = open(output_file, "w")
        GFF.write(GFF.parse(in_fd, target_lines=target_lines), out_fd)
        in_fd.close()
        out_fd.close()

    @staticmethod
    def gff32gtf(input_file, output_file):
        os.system("gffread %s -T -o %s" % (input_file, output_file))

    @staticmethod
    def gtf2gff3(input_file, output_file):
        os.system("gffread %s -o %s" % (input_file, output_file))


