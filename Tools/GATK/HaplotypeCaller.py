#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool

import os
from Routines.Functions import check_path


class HaplotypeCaller(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T HaplotypeCaller", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    @staticmethod
    def parse_options_for_parallel_run(reference, alignment, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                                       stand_call_conf=30, gvcf_mode=False):

        options = " -R %s" % reference
        options += " -I %s" % alignment
        options += " --genotyping_mode %s" % genotyping_mode if genotyping_mode else ""
        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " -stand_call_conf %i" % stand_call_conf
        options += " --emitRefConfidence GVCF" if gvcf_mode else ""

        return options

    def parse_options(self, reference, alignment, output, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                      stand_call_conf=30, gvcf_mode=False, include_region_id_file=None, exclude_region_id_file=None):

        options = " -nct %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % alignment
        options += " --genotyping_mode %s" % genotyping_mode if genotyping_mode else ""
        options += " --output_mode %s" % output_mode if output_mode else ""
        #options += " -stand_emit_conf %i" % stand_emit_conf
        options += " -stand_call_conf %i" % stand_call_conf
        options += " --emitRefConfidence GVCF" if gvcf_mode else ""
        options += " -L %s" % include_region_id_file if include_region_id_file else ""
        options += " -XL %s" % exclude_region_id_file if exclude_region_id_file else ""

        options += " -o %s" % output

        return options

    def call(self, reference, alignment, output, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
             stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
            java -Xmx100g -jar ~/tools/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
              -T HaplotypeCaller \
              -R ${fasta} \
              -I ${bam%bam}realigned.bam \
              --genotyping_mode DISCOVERY \
              --output_mode EMIT_VARIANTS_ONLY \
              -stand_call_conf 30 \
              -o ${bam%bam}raw.vcf
        """
        options = self.parse_options(reference, alignment, output, genotyping_mode=genotyping_mode,
                                     output_mode=output_mode,
                                     stand_call_conf=stand_call_conf, gvcf_mode=False,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options)

    def gvcf_call(self, reference, alignment, output, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                  stand_call_conf=30, include_region_id_file=None, exclude_region_id_file=None):
        """
        java -jar GenomeAnalysisTK.jar \
         -R reference.fasta \
         -T HaplotypeCaller \
         -I sample1.bam \
         --emitRefConfidence GVCF \
         [--dbsnp dbSNP.vcf] \
         [-L targets.interval_list] \
         -o output.raw.snps.indels.g.vcf
        """
        options = self.parse_options(reference, alignment, output, genotyping_mode=genotyping_mode,
                                     output_mode=output_mode,
                                     stand_call_conf=stand_call_conf, gvcf_mode=True,
                                     include_region_id_file=include_region_id_file,
                                     exclude_region_id_file=exclude_region_id_file)

        self.execute(options)

    def parallel_gvcf_call(self, reference, alignment, output_dir, output_prefix,
                           genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                           stand_call_conf=30, max_region_length=1000000, max_seqs_per_region=500,
                           length_dict=None, parsing_mode="parse", region_list=None, ):
        self.safe_mkdir(output_dir)

        region_list = self.prepare_region_list_by_length(max_length=max_region_length,
                                                         max_seq_number=max_seqs_per_region,
                                                         length_dict=length_dict,
                                                         reference=None if length_dict is not None else reference,
                                                         parsing_mode=parsing_mode,
                                                         output_dir="%s/regions/" % output_dir) if region_list is None else region_list

        options = self.parse_options_for_parallel_run(reference, alignment,
                                                      genotyping_mode=genotyping_mode,
                                                      output_mode=output_mode,
                                                      stand_call_conf=stand_call_conf,
                                                      gvcf_mode=True)
        options += " -nct 1"
        options_list = []

        output_index = 1
        for regions in region_list:

            region_options = " -o %s/%s_%i.g.vcf" % (output_dir, output_prefix, output_index)
            for region in regions:
                region_options += " -L %s:%i-%i" % (region[0], region[1], region[2])

            options_list.append(options + region_options)
            output_index += 1

        self.parallel_execute(options_list)

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
        """
        OBSOLETE FUNCTION!!! DOESN'T WORK
        """
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
        os.system(" java -jar %sGenomeAnalysisTK.jar -nt %i -l INFO -R %s -T HaplotypeCaller -I %s -stand_call_conf %i -stand_emit_conf %i -o %s --output_mode %s -glm %s %s"
                  % (gatk_dir, num_of_threads, reference_file, alignment, stand_call_conf, stand_emit_conf, output_file, output_mode, discovery_mode, default_qualities))