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
                     output_file="GATK_raw.vcf",
                     default_base_qualities=None):
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
        default_qualities = ""
        if default_base_qualities:
            default_qualities = "--defaultBaseQualities %i" % default_base_qualities

        gatk_dir = check_path(GATK_dir)
        os.system(" java -jar %sGenomeAnalysisTK.jar -nt %i -l INFO -R %s -T UnifiedGenotyper -I %s -stand_call_conf %i -stand_emit_conf %i -o %s --output_mode %s -glm %s %s"
                  % (gatk_dir, num_of_threads, reference_file, alignment, stand_call_conf, stand_emit_conf, output_file, output_mode, discovery_mode, default_qualities))


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
    # TODO: fix problems with unrecognized data
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.html

    def build_model(self, reference_file, input_vcf, training_set_good,
                       training_set_bad=None, mode="SNP",
                       recalfile="output_recal",
                       tranchesfile="output.tranches",
                       #rscriptfile="output.plots.R",
                       GATK_dir=""):
        # possible mode:
        # SNP
        # INDEL
        bad_set = ""
        if training_set_bad:
            bad_set = "-resource:bad,known=false,training=true,bad=true %s" % training_set_bad
            """
            print("java -Xmx4g -jar %sGenomeAnalysisTK.jar -T VariantRecalibrator -R %s -input %s -an QD \
                  -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
                  -mode %s -resource:good,training=true,truth=true %s %s -recalFile %s -tranchesFile %s -rscriptFile %s"
                  % (GATK_dir, reference_file, input_vcf, mode, training_set_good, bad_set, recalfile, tranchesfile, rscriptfile))
            """
        os.system("java -Xmx4g -jar %sGenomeAnalysisTK.jar -T VariantRecalibrator -R %s -input %s -an QD \
                  -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
                  -mode %s -resource:good,known=true,training=true,truth=true %s %s -recalFile %s -tranchesFile %s"
                  % (GATK_dir, reference_file, input_vcf, mode, training_set_good, bad_set, recalfile, tranchesfile))


class ApplyRecalibration():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_ApplyRecalibration.html

    def apply_model(self, reference_file, input_vcf, tranchesfile, recalfile,
                    outputfile, mode="SNP", ts_filter_level=99.0, GATK_dir=""):
        # possible mode:
        # SNP
        # INDEL
        os.system("java -Xmx3g -jar %sGenomeAnalysisTK.jar -T ApplyRecalibration -R %s -input %s \
                   --ts_filter_level %f -tranchesFile %s -recalFile %s -mode %s -o %s"
                  % (GATK_dir, reference_file, input_vcf, ts_filter_level, tranchesfile, recalfile, mode, outputfile))


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

    def combine_from_different_sources(self, gatk_dir, reference, output_file,
                                       file_list, merge_option="REQUIRE_UNIQUE"):
        variants_string = " ".join(list(map(lambda x: "--variant %s" % x, file_list)))

        os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T CombineVariants %s -o %s -genotypeMergeOptions %s"
                  % (gatk_dir, reference, variants_string, output_file, merge_option))

if __name__ == "__main__":
    gatk_dir = "/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0/"
    reference_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.fasta"
    """
    input_vcf = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/N012-LAN210-Can-PmCDA1-NA-RUN2-D3/alignment_LAN210_v0.9m/N012-LAN210-Can-PmCDA1-NA-RUN2-D3_GATK_raw.vcf"
    training_set_good = "/home/mahajrod/genetics/desaminases/data/model_vcf/LEFT/210-L1_GATK_best_snps.recode_LEFT.vcf"
    training_set_bad = "/home/mahajrod/genetics/desaminases/data/model_vcf/BAD/210-L1_GATK_best_snps.recode_BAD.vcf"
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/N011-LAN210-Can-PmCDA1-NA-RUN2-D3/alignment_LAN210_v0.9m/recalibrate/"
    recalfile = "%soutput_recal" % workdir
    tranchesfile = "%soutput.tranches" % workdir
    rscriptfile = "%soutput.plots.R" % workdir
    outputfile = "%recalibrated.vcf" % workdir

    var_recalibrator = VariantRecalibrator()
    app_recalibration = ApplyRecalibration()

    var_recalibrator.build_model(reference_file, input_vcf, training_set_good,
                               training_set_bad=training_set_bad, mode="SNP",
                               recalfile=recalfile,
                               tranchesfile=tranchesfile,
                               #rscriptfile=rscriptfile,
                               GATK_dir=gatk_dir)

    app_recalibration = ApplyRecalibration()

    app_recalibration.apply_model(reference_file, input_vcf, tranchesfile, recalfile,
                    outputfile, mode="SNP", ts_filter_level=99.0, GATK_dir=gatk_dir)
    """
    file_list = ["/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/N010-LAN210-Can-PmCDA1-NA-RUN2-D3/alignment_LAN210_v0.9m/N010-LAN210-Can-PmCDA1-NA-RUN2-D3_GATK_best_merged.vcf",
                 "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/N011-LAN210-Can-PmCDA1-NA-RUN2-D3/alignment_LAN210_v0.9m/N011-LAN210-Can-PmCDA1-NA-RUN2-D3_GATK_best_merged.vcf"]
    combo = CombineVariants()

    combo.combine_from_different_sources(gatk_dir, reference_file, "out_file.vcf", file_list)

