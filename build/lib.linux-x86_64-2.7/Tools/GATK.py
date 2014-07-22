#!/usr/bin/env python
import os
from General.General import check_path


class UnifiedGenotyper():

    def variant_call(self,
                     alignment,
                     reference_file,
                     stand_emit_conf=40,
                     stand_call_conf=100,
                     GATK_dir="",
                     num_of_threads=5,
                     output_mode="EMIT_VARIANTS_ONLY",
                     discovery_mode="BOTH",
                     output_file="GATK_raw.vcf"):
        # manual http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.html
        # output_mode values:
        # EMIT_VARIANTS_ONLY
        # EMIT_ALL_CONFIDENT_SITES
        # EMIT_ALL_SITES

        # discovery_mode values:
        # SNP
        # INDEL
        # GENERALPLOIDYSNP
        # GENERALPLOIDYINDEL
        # BOTH

        gatk_dir = check_path(GATK_dir)
        os.system(" java -jar %sGenomeAnalysisTK.jar -nt %i -l INFO -R %s -T UnifiedGenotyper -I %s -stand_call_conf %i -stand_emit_conf %i -o %s --output_mode %s -glm %s"
                  % (gatk_dir, num_of_threads, reference_file, alignment, stand_call_conf, stand_emit_conf, output_file, output_mode, discovery_mode))


class SelectVariants():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_SelectVariants.html
    # selectType
    # INDEL
    # SNP
    # MIXED
    # MNP
    # SYMBOLIC
    # NO_VARIATION
    def select_variants(self, gatk_dir, reference_file, input_vcf, output_vcf, vartype=None, varfilter=None):
        selecttype = ""
        if vartype:
            selecttype = "-selectType \'%s\'" % vartype

        filter_exp = ""
        if varfilter:
            filter_exp = "-select \'%s\'" % varfilter
        os.system("java -jar %sGenomeAnalysisTK.jar -T SelectVariants -R %s -V %s %s %s -o %s"
                  % (gatk_dir, reference_file, input_vcf, selecttype, filter_exp, output_vcf))

    def get_SNP(self, gatk_dir, reference_file, input_vcf, output_vcf):
        self.select_variants(gatk_dir, reference_file, input_vcf, output_vcf, vartype="SNP")

    def get_indel(self, gatk_dir, reference_file, input_vcf, output_vcf):
        self.select_variants(gatk_dir, reference_file, input_vcf, output_vcf, vartype="INDEL")

    def remove_filtered(self, gatk_dir, reference_file, input_vcf, output_vcf):
        os.system("java -jar %sGenomeAnalysisTK.jar -T SelectVariants -R %s -V %s -o %s -ef"
                  % (gatk_dir, reference_file, input_vcf, output_vcf))

    #Use STR filtration based on GATK prediction very carefully as it counts even AAA -> AA as STR
    def get_STR(self, gatk_dir, reference_file, input_vcf, output_vcf):
        self.select_variants(gatk_dir, reference_file, input_vcf, output_vcf, varfilter="STR")

    def get_nonSTR(self, gatk_dir, reference_file, input_vcf, output_vcf):
        self.select_variants(gatk_dir, reference_file, input_vcf, output_vcf, varfilter="vc.hasAttribute(\"STR\") == 0")


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


class VariantRecalibrator():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.html

    def filter_bad_SNP(self):
        pass


class FastaAlternateReferenceMaker():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html

    def correct_reference(self, gatk_dir, reference, new_reference, variants_vcf):

        os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --variant %s"
                  % (gatk_dir, reference, new_reference, variants_vcf))


class CombineVariants():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html

    def combine_from_same_source(self, gatk_dir, reference, input_filelist, output_file, num_threads=4):
        in_files = ""
        for entry in input_filelist:
            in_files += " --variant %s" % entry
        os.system("java -jar %sGenomeAnalysisTK.jar -R %s -T CombineVariants -nt %i %s -o %s"
                  % (gatk_dir, reference, num_threads, in_files, output_file))
