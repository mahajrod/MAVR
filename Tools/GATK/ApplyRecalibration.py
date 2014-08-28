#!/usr/bin/env python
import os
from Tools.Abstract import Tool


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
