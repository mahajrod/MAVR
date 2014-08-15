import os
from Routines.Functions import check_path


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