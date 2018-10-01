import os
from Routines.Functions import check_path
from Tools.Abstract import JavaTool


class UnifiedGenotyper(JavaTool):

    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T UnifiedGenotyper", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def variant_call(self,
                     alignment,
                     reference_file,
                     stand_emit_conf=40,
                     stand_call_conf=100,
                     output_mode="EMIT_VARIANTS_ONLY",
                     discovery_mode="BOTH",
                     output_file="GATK_raw.vcf",
                     default_base_qualities=None,
                     quality_to_change=None,
                     quality_to_change_to=60,
                     remove_spliced_reads=None,
                     use_spliced_reads=None):


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

        options = " -nt %i" % self.threads
        options += " -R %s" % reference_file

        options += " -I %s" % alignment
        options += " -stand_call_conf %i" % stand_call_conf
        options += " -stand_emit_conf %i" % stand_emit_conf
        options += " -o %s" % output_file
        options += " --output_mode %s" % output_mode
        options += " -glm %s" % discovery_mode
        options += " --defaultBaseQualities %i" % default_base_qualities if default_base_qualities else ""
        options += " --filter_reads_with_N_cigar" if remove_spliced_reads else ""
        options += " -U ALLOW_N_CIGAR_READS" if use_spliced_reads else ""

        if quality_to_change and quality_to_change_to:
            options += " -rf ReassignOneMappingQuality -RMQF %i -RMQT %i" % (quality_to_change, quality_to_change_to)

        self.execute(options=options)

        """
        default_qualities = ""
        if default_base_qualities:
            default_qualities = "--defaultBaseQualities %i" % default_base_qualities

        #gatk_dir = check_path(GATK_dir)
        os.system(" java -jar %sGenomeAnalysisTK.jar -nt %i -l INFO -R %s -T UnifiedGenotyper -I %s -stand_call_conf %i -stand_emit_conf %i -o %s --output_mode %s -glm %s %s"
                  % (gatk_dir, num_of_threads, reference_file, alignment, stand_call_conf, stand_emit_conf, output_file, output_mode, discovery_mode, default_qualities))
        """
