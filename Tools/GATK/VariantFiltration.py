#!/usr/bin/env python
import os
from Tools.Abstract import JavaTool



class VariantFiltration(JavaTool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html
    # default filters for indel and snp filtration were taken from GATK BestPractice
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
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
                       MQ=40.0, HaplotypeScore=13.0, MappingQualityRankSum=-12.5, ReadPosRankSum=-8.0):
        """
        filter_expression = "'QD < %f || FS > %f || MQ < %f || HaplotypeScore > %f || MQRankSum < %f || ReadPosRankSum < %f'" \
                           % (QD, FS, MQ, HaplotypeScore, MappingQualityRankSum, ReadPosRankSum)
        """
        filter_expression = "'QD < %f || FS > %f || MQ < %f || MQRankSum < %f || ReadPosRankSum < %f'" \
                           % (QD, FS, MQ, MappingQualityRankSum,  ReadPosRankSum)

        self.filter(reference_file, input_vcf, output_vcf, filter_expression, filter_name)

    def filter_bad_indel(self, reference_file, input_vcf, output_vcf, filter_name='ambiguous_indel', QD=2.0,
                         ReadPosRankSum=-20.0, InbreedingCoeff=-0.8, FS=200.0):
        """
        filter_expression = "'QD < %f || ReadPosRankSum < %f || InbreedingCoeff < %f || FS > %f'" \
                           % (QD, ReadPosRankSum, InbreedingCoeff, FS)
        """
        filter_expression = "'QD < %f || ReadPosRankSum < %f || FS > %f'" \
                           % (QD, ReadPosRankSum, FS)
        self.filter(reference_file, input_vcf, output_vcf, filter_expression, filter_name)

    def filter_bad_variants(self, reference_file, input_vcf, output_prefix, snp_filter_name='ambiguous_snp', snp_QD=2.0,
                            snp_FS=60.0, snp_MQ=40.0, snp_HaplotypeScore=13.0, snp_MappingQualityRankSum=-12.5,
                            snp_ReadPosRankSum=-8.0, indel_filter_name='ambiguous_indel', indel_QD=2.0,
                            indel_ReadPosRankSum=-20.0, indel_InbreedingCoeff=-0.8, indel_FS=200.0):
        from Tools.GATK import SelectVariants, CombineVariants
        snp_raw_vcf = "%s.snp.raw.vcf" % output_prefix
        indel_raw_vcf = "%s.indel.raw.vcf" % output_prefix

        snp_filtered_vcf = "%s.snp.filtered.vcf" % output_prefix
        indel_filtered_vcf = "%s.indel.filtered.vcf" % output_prefix

        combined_filtered_vcf = "%s.combined.filtered.vcf" % output_prefix
        print type(SelectVariants)
        SelectVariants.jar_path = self.jar_path
        CombineVariants.jar_path = self.jar_path
        SelectVariants.get_SNP(reference_file, input_vcf, snp_raw_vcf)
        SelectVariants.get_indel(reference_file, input_vcf, indel_raw_vcf)

        self.filter_bad_SNP(reference_file, snp_raw_vcf, snp_filtered_vcf, filter_name=snp_filter_name, QD=snp_QD,
                            FS=snp_FS, MQ=snp_MQ, HaplotypeScore=snp_HaplotypeScore,
                            MappingQualityRankSum=snp_MappingQualityRankSum, ReadPosRankSum=snp_ReadPosRankSum)
        self.filter_bad_indel(reference_file, indel_raw_vcf, indel_filtered_vcf, filter_name=indel_filter_name,
                              QD=indel_QD, ReadPosRankSum=indel_ReadPosRankSum, InbreedingCoeff=indel_InbreedingCoeff,
                              FS=indel_FS)
        CombineVariants.combine_from_same_source(reference_file, [snp_filtered_vcf, indel_filtered_vcf],
                                                 combined_filtered_vcf)

    """
    def filter(self, gatk_dir, reference_file, input_vcf, filter_expresion, filter_name, output_vcf):
        #print("java -jar %sGenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s --filterExpression \'%s\' --filterName %s -o %s"
        #          % (gatk_dir, reference_file, input_vcf, filter_expresion, filter_name, output_vcf))
        os.system("java -jar %sGenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s --filterExpression \'%s\' --filterName %s -o %s"
                  % (gatk_dir, reference_file, input_vcf, filter_expresion, filter_name, output_vcf))

    def filter_bad_SNP(self, gatk_dir, reference_file, input_vcf, output_vcf, QD=2.0, FS=60.0, MQ=40.0,
                       HaplotypeScore=13.0, MappingQualityRankSum=-12.5, ReadPosRankSum=-8.0):
        filter_expresion = 'QD < %f || FS > %f || MQ < %f || HaplotypeScore > %f || MappingQualityRankSum < %f || ReadPosRankSum < %f' \
                           % (QD, FS, MQ, HaplotypeScore, MappingQualityRankSum, ReadPosRankSum)
        print(filter_expresion)
        filter_name = 'ambigious_snp'
        self.filter(gatk_dir, reference_file, input_vcf, filter_expresion, filter_name, output_vcf)

    def filter_bad_indel(self, gatk_dir, reference_file, input_vcf, output_vcf, QD=2.0,
                         ReadPosRankSum=-20.0, InbreedingCoeff=-0.8, FS=200.0):
        filter_expresion = "QD < %f || ReadPosRankSum < %f || InbreedingCoeff < %f || FS > %f" \
                           % (QD, ReadPosRankSum, InbreedingCoeff, FS)
        filter_name = 'ambigious_indel'
        self.filter(gatk_dir, reference_file, input_vcf, filter_expresion, filter_name, output_vcf)
    """

