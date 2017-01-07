#!/usr/bin/env python
import os
from Tools.Alignment import Bowtie2

from Pipelines.Abstract import Pipeline


class ScaffoldingPipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self)

    def prepare_directories_for_insert_size_estimation(self, output_directory, sample_list):
        out_dir = os.path.abspath(output_directory)
        index_directory = "%s/bowtie_index/" % out_dir

        self.save_mkdir(index_directory)
        for sample in sample_list:
            sample_dir = "%s/%s" % (out_dir, sample)
            self.save_mkdir(sample_dir)

    def get_insert_size_distribution(self, sample_directory, forward_files, reverse_files, estimated_insert_size,
                                     output_prefix, genome, genome_index, read_orientation="fr",
                                     parsing_mode="index_db", number_of_bins=100, genome_format="fasta"):

        sample_dir = os.path.abspath(sample_directory)
        output_pref = "%s/%s" % (sample_dir, output_prefix)
        min_contig_len_threshold = 3 * estimated_insert_size
        region_bed_file = "%s/contig.bed" % sample_dir
        self.make_region_bed_file_from_file(genome, region_bed_file, min_len=min_contig_len_threshold,
                                            parsing_mode=parsing_mode, input_format=genome_format)

        output_len_file = "%s.len" % output_pref

        bowtie_string = "bowtie2 --very-sensitive -x %s -1 %s -2 %s  -p %i -X %i --%s " % (genome_index,
                                                                                           forward_files if isinstance(forward_files, str) else ",".join(forward_files),
                                                                                           reverse_files if isinstance(reverse_files, str) else ",".join(reverse_files),
                                                                                           self.threads,
                                                                                           min_contig_len_threshold,
                                                                                           read_orientation)

        samtools_string = "samtools view -L %s -" % region_bed_file
        awk_string = "awk -F'\\t' '{ if ($9 > 0) print $9}'"

        final_string = "%s | %s | %s > %s" % (bowtie_string, samtools_string, awk_string, output_len_file)

        self.execute(cmd=final_string)

        self.draw_histogram(output_len_file, output_pref, max_length=min_contig_len_threshold,
                            number_of_bins=number_of_bins, xlabel="Insert size",
                            ylabel="Number of fragments", title="Insert size distribution", extensions=("png",))
