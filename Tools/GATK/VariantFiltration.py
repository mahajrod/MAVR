#!/usr/bin/env python
import os


class VariantFiltration():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html
    # default filters for indel and snp filtration were taken from GATK BestPractice
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

