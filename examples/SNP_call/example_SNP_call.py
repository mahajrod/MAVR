#!/usr/bin/env python
from Pipelines.SNPCall import *
from Tools.Picard import *
from Tools.GATK import UnifiedGenotyper, SelectVariants, VariantFiltration
#in case reads are already trimmed and corrected
reference = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.3m/LAN210_v0.3m_masked.fasta"
reference_dict = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.3m/LAN210_v0.3m_masked.dict"
reference_index = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.3m/LAN210_v0.3m_masked"
gatk_dir = check_path("/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.1-1")
run_name = "Original_reads"
work_dir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.7m_masked_original"

sample_name = "N008-LAN210-wt-NA-NA-RUN2"
forward_reads = "/home/mahajrod/genetics/desaminases/data/LAN210_raw/LAN210_corrected/210_ACTTGA_L003_edited_R_1.fastq"
reverse_reads = "/home/mahajrod/genetics/desaminases/data/LAN210_raw/LAN210_corrected/210_ACTTGA_L003_edited_R_2.fastq"
platform = "illumina"
os.chdir(work_dir)
#os.system("source /etc/profile")
"""
get_chromosomes_bed(reference, reference_index, mitochondrial_region_name="mt",
                        chrom_out_file="chromosomes.bed", mito_out_file="mt.bed", reference_filetype="fasta")

get_alignment_without_trim(reference_index,
                              sample_name,
                              forward_reads,
                              "chromosomes.bed",
                              "mt.bed",
                              reverse_reads=reverse_reads,
                              max_threads=5,
                              quality_score="phred33")

add_header2bam("%s_trimmed_sorted_rm_pcr_chrom.bam" % sample_name,
                "%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                sample_name,
                sample_name,
                platform,
                sample_name,
                sample_name,
                PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115")

make_fasta_dict(reference, reference_dict, PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115")

os.system("samtools faidx %s" % reference)
snp_call_GATK("%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                sample_name,
                reference,
                "", #ignored if skip base recall is True
                stand_call_conf=100,
                GATK_dir="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.1-1",
                skip_base_recalibration=True)

#remove filtered SNPs -ef flag

os.system("java -jar %sGenomeAnalysisTK.jar -R %s -T SelectVariants --variant %s_GATK_filtered_snps.vcf -o %s_GATK_best_snps.vcf -ef"
          % (gatk_dir, reference, sample_name, sample_name))
"""
vcf_all = "%s_GATK_raw.vcf" % sample_name
vcf_indel = "%s_GATK_raw_indel.vcf" % sample_name
vcf_SNP = "%s_GATK_raw_SNP.vcf" % sample_name
vcf_filtered_indel = "%s_GATK_filtered_indel.vcf" % sample_name
vcf_filtered_SNP = "%s_GATK_filtered_SNP.vcf" % sample_name
vcf_best_indel = "%s_GATK_best_indel.vcf" % sample_name
vcf_best_SNP = "%s_GATK_best_SNP.vcf" % sample_name

UnifiedGenotyper.variant_call("%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                              reference,
                              stand_emit_conf=40,
                              stand_call_conf=100,
                              GATK_dir=gatk_dir,
                              num_of_threads=5,
                              output_mode="EMIT_VARIANTS_ONLY",
                              discovery_mode="BOTH",
                              output_file=vcf_all)

SelectVariants.get_indel(gatk_dir, reference, vcf_all, vcf_indel)
SelectVariants.get_SNP(gatk_dir, reference, vcf_all, vcf_SNP)

VariantFiltration.filter_bad_SNP(gatk_dir, reference, vcf_SNP, vcf_filtered_SNP)
VariantFiltration.filter_indel(gatk_dir, reference, vcf_indel, vcf_filtered_indel)

VariantFiltration.remove_filtered(gatk_dir, reference, vcf_filtered_SNP, vcf_best_SNP)
VariantFiltration.remove_filtered(gatk_dir, reference, vcf_filtered_indel, vcf_best_indel)

#os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --variant %s_GATK_best_snps.vcf"
#          % (gatk_dir, reference, "LAN210_v0.7m.fasta", sample_name))

#os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --snpmask mask.vcf %s_GATK_best_snps.vcf"
#          % (gatk_dir, reference, "LAN210_v0.7m_masked_SNP.fasta", sample_name))

