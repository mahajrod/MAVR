#!/usr/bin/env python
import os


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

