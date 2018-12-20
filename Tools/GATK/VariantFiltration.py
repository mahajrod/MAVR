#!/usr/bin/env python

from Tools.Abstract import JavaTool
from Routines import VCFRoutines


class VariantFiltration(JavaTool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html
    # default filters for indel and snp filtration were taken from GATK BestPractice
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T VariantFiltration", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    @staticmethod
    def parse_common_options(reference_file, input_vcf, output_vcf, filter_expression, filter_name):
        options = " -R %s" % reference_file
        options += " -V %s" % input_vcf
        options += " --filterExpression \'%s\'" % filter_expression
        options += " --filterName %s" % filter_name
        options += " -o %s" % output_vcf

        return options

    def filter(self, reference_file, input_vcf, output_vcf, filter_expression, filter_name, ):
        options = self.parse_common_options(reference_file, input_vcf, output_vcf, filter_expression, filter_name)

        self.execute(options=options)

    def filter_bad_SNP(self, reference_file, input_vcf, output_vcf, filter_name='ambiguous_snp', QD=2.0, FS=60.0,
                       MQ=40.0, MappingQualityRankSum=-12.5, ReadPosRankSum=-8.0):
        """
        "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        """
        filter_expression = "QD < %f || FS > %f || MQ < %f || MQRankSum < %f || ReadPosRankSum < %f" \
                           % (QD, FS, MQ, MappingQualityRankSum,  ReadPosRankSum)

        self.filter(reference_file, input_vcf, output_vcf, filter_expression, filter_name)

    def filter_bad_indel(self, reference_file, input_vcf, output_vcf, filter_name='ambiguous_indel', QD=2.0,
                         ReadPosRankSum=-20.0,  FS=200.0):
        """
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        """
        filter_expression = "QD < %f || ReadPosRankSum < %f || FS > %f" % (QD, ReadPosRankSum, FS)
        self.filter(reference_file, input_vcf, output_vcf, filter_expression, filter_name)

    def filter_bad_variants(self, reference_file, input_vcf, output_prefix, snp_filter_name='ambiguous_snp', snp_QD=2.0,
                            snp_FS=60.0, snp_MQ=40.0, snp_HaplotypeScore=13.0, snp_MappingQualityRankSum=-12.5,
                            snp_ReadPosRankSum=-8.0, indel_filter_name='ambiguous_indel', indel_QD=2.0,
                            indel_ReadPosRankSum=-20.0, indel_FS=200.0, combine_vcf=False, sequence_dict_file=None,
                            picard_memory="1g", picard_dir=None):

        from Tools.GATK import SelectVariants
        from Tools.Picard import SortVcf
        snp_raw_vcf = "%s.snp.raw.vcf" % output_prefix
        indel_raw_vcf = "%s.indel.raw.vcf" % output_prefix

        snp_filtered_vcf = "%s.snp.with_filters.vcf" % output_prefix
        indel_filtered_vcf = "%s.indel.with_filters.vcf" % output_prefix

        snp_good_vcf = "%s.snp.good.vcf" % output_prefix
        indel_good_vcf = "%s.indel.good.vcf" % output_prefix

        unsorted_combined_filtered_vcf = "%s.combined.with_filters.unsorted.vcf" % output_prefix
        unsorted_combined_good_vcf = "%s.combined.good.unsorted.vcf" % output_prefix

        combined_filtered_vcf = "%s.combined.with_filters.sorted.vcf" % output_prefix
        combined_good_vcf = "%s.combined.good.sorted.vcf" % output_prefix


        SelectVariants.jar_path = self.jar_path
        SortVcf.jar_path = picard_dir
        SortVcf.max_memory = picard_memory
        #CombineVariants.jar_path = self.jar_path

        SelectVariants.get_SNP(reference_file, input_vcf, snp_raw_vcf)
        SelectVariants.get_indel(reference_file, input_vcf, indel_raw_vcf)

        self.filter_bad_SNP(reference_file, snp_raw_vcf, snp_filtered_vcf, filter_name=snp_filter_name, QD=snp_QD,
                            FS=snp_FS, MQ=snp_MQ, #HaplotypeScore=snp_HaplotypeScore,
                            MappingQualityRankSum=snp_MappingQualityRankSum, ReadPosRankSum=snp_ReadPosRankSum)
        self.filter_bad_indel(reference_file, indel_raw_vcf, indel_filtered_vcf, filter_name=indel_filter_name,
                              QD=indel_QD, ReadPosRankSum=indel_ReadPosRankSum, FS=indel_FS)

        SelectVariants.remove_entries_with_filters(reference_file, snp_filtered_vcf, snp_good_vcf)
        SelectVariants.remove_entries_with_filters(reference_file, indel_filtered_vcf, indel_good_vcf)

        if combine_vcf:
            VCFRoutines.combine_same_samples_vcfs(unsorted_combined_filtered_vcf, vcf_list=[snp_filtered_vcf, indel_filtered_vcf],
                                                  order_vcf_files=False, sort=True, chunk_folder=None, chunk_prefix=None,
                                                  chunk_suffix=None, starting_chunk=None, chunk_number_list=None,
                                                  close_fd_after=False, extension_list=[".vcf", ])

            VCFRoutines.combine_same_samples_vcfs(unsorted_combined_good_vcf, vcf_list=[snp_good_vcf, indel_good_vcf],
                                                  order_vcf_files=False, sort=True, chunk_folder=None, chunk_prefix=None,
                                                  chunk_suffix=None, starting_chunk=None, chunk_number_list=None,
                                                  close_fd_after=False, extension_list=[".vcf", ])
            if sequence_dict_file:
                SortVcf.sort_vcf(unsorted_combined_filtered_vcf, combined_filtered_vcf, seq_dict=sequence_dict_file)
                SortVcf.sort_vcf(unsorted_combined_good_vcf, combined_good_vcf, seq_dict=sequence_dict_file)
            """
            #CombineVariants IS TOO SLOW!!!! It takes DAYS to merge VCFs
            CombineVariants.combine_from_same_source(reference_file, [snp_filtered_vcf, indel_filtered_vcf],
                                                     combined_filtered_vcf)

            CombineVariants.combine_from_same_source(reference_file, [snp_good_vcf, indel_good_vcf],
                                                     combined_good_vcf)
            """
