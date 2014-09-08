import os
from Tools.Abstract import Tool


class SelectVariants(Tool):
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

