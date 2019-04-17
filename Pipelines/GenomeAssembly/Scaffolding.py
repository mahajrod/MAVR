#!/usr/bin/env python
import os
from RouToolPa.Tools.Samtools import SamtoolsV1
from Pipelines.Abstract import Pipeline



class ScaffoldingPipeline(Pipeline):
    def __init__(self):
        Pipeline.__init__(self)

    def prepare_directories_for_insert_size_estimation(self, output_directory, sample_list):
        out_dir = os.path.abspath(output_directory)
        index_directory = "%s/bowtie_index/" % out_dir

        self.safe_mkdir(index_directory)
        for sample in sample_list:
            sample_dir = "%s/%s" % (out_dir, sample)
            self.safe_mkdir(sample_dir)

    def get_insert_size_distribution(self, sample_directory, forward_files, reverse_files, estimated_insert_size,
                                     output_prefix, genome, genome_index, input_files_are_fasta=False,
                                     read_orientation="fr", parsing_mode="index_db", number_of_bins=100,
                                     genome_format="fasta", store_sam=False, aligner="bowtie2",
                                     aligner_binary_dir="", xlimit_for_histo=None):

        sample_dir = os.path.abspath(sample_directory)
        output_pref = "%s/%s" % (sample_dir, output_prefix)
        min_contig_len_threshold = 3 * estimated_insert_size
        region_bed_file = "%s/%s.contig.bed" % (sample_dir, output_prefix)
        self.make_region_bed_file_from_file(genome, region_bed_file, min_len=min_contig_len_threshold,
                                            parsing_mode=parsing_mode, input_format=genome_format)

        output_filtered_len_file = "%s.filtered.len" % output_pref
        output_all_len_file = "%s.all.len" % output_pref
        output_sam = "%s.sam" % output_pref
        output_bam = "%s.bam" % output_pref
        aligner_log = "%s.aligner.log"
        forward_reads = forward_files if isinstance(forward_files, str) else ",".join(forward_files)
        reverse_reads = reverse_files if isinstance(reverse_files, str) else ",".join(reverse_files)

        bowtie_options = " --very-sensitive"
        bowtie_options += " -x %s" % genome_index
        bowtie_options += " -1 %s" % forward_reads

        bowtie_options += " -2 %s" % reverse_reads
        bowtie_options += " -p %i" % self.threads
        bowtie_options += " -X %i" % min_contig_len_threshold
        bowtie_options += " --%s" % read_orientation
        bowtie_options += " -f" if input_files_are_fasta else ""

        bowtie2_string = "bowtie2 %s 2>%s" % (bowtie_options, aligner_log)

        bwa_options = " mem"
        bwa_options += " -t %i" % self.threads
        bwa_options += " %s" % genome_index
        bwa_options += " %s %s" % (forward_reads, reverse_reads)

        bwa_string = "bwa %s" % bwa_options

        tee_string = "tee %s" % output_sam
        samtools_string = "samtools view -L %s -" % region_bed_file
        awk_string = "awk -F'\\t' '{ if ($9 > 0) print $9}'"

        if aligner == "bowtie2":
            aligner_string = bowtie2_string
        elif aligner == "bwa":
            aligner_string = bwa_string
        else:
            raise ValueError("Unrecognized aligner: %s" % aligner)

        aligner_string = "%s%s" % (self.check_path(aligner_binary_dir), aligner_string)

        final_string = "%s | %s | %s | %s > %s" % (aligner_string, tee_string, samtools_string,
                                                   awk_string, output_filtered_len_file)

        full_len_string = "%s %s > %s" % (awk_string, output_sam, output_all_len_file)

        self.execute(cmd=final_string)
        self.execute(cmd=full_len_string)

        self.draw_histogram_from_file(output_filtered_len_file, output_pref,
                            max_length=xlimit_for_histo if xlimit_for_histo else min_contig_len_threshold,
                            number_of_bins=number_of_bins, xlabel="Insert size",
                            ylabel="Number of fragments", title="Insert size distribution", extensions=("png", "svg"))

        SamtoolsV1.convert_sam_and_index(output_sam, output_bam)

        if not store_sam:
            os.remove(output_sam)


    #def get
