#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class TVC(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "tvc", path=path, max_threads=max_threads)

    def variant_caller_pipeline(self, input_bam, reference_fasta, out_dir="./", hotspot_vcf=None, param_file=None,
                                tvc_root_dir=None, primer_trim_bed=None, postprocesed_bam=None, region_bed=None,
                                error_motif_file=None):
        options = ""
        options += " -b %s" % region_bed if region_bed else ""
        options += " -s %s" % hotspot_vcf if hotspot_vcf else ""
        options += " -i %s" % input_bam
        options += " -r %s" % reference_fasta
        options += " -o %s" % out_dir
        options += " -p %s" % param_file if param_file else ""
        options += " -B %s" % tvc_root_dir if tvc_root_dir else ""
        options += " -n %i" % self.threads
        options += " --primer-trim-bed=%s" % primer_trim_bed if primer_trim_bed else ""
        options += " --postprocessed-bam=%s" % postprocesed_bam if postprocesed_bam else ""


        self.execute(options, cmd="variant_caller_pipeline.py")

    def variant_call(self, input_bam, reference_fasta, out_vcf, out_dir="./", sample_name=None, num_variants_per_thread=None,
                     parameters_file=None, force_sample_name=None, target_file=None, trim_ampliseq_primers=None,
                     downsample_to_coverage=None, model_file=None, recal_model_hp_thres=None,
                     output_vcf=None, suppress_reference_genotypes=None, suppress_no_calls=None,
                     suppress_nocall_genotypes=None, allow_snps=None, allow_indels=None, allow_mnps=None,
                     allow_complex=None, max_complex_gap=None, use_best_n_alleles=None, min_mapping_qv=None,
                     read_snp_limit=None, read_max_mismatch_fraction=None, gen_min_alt_allele_freq=None,
                     gen_min_indel_alt_allele_freq=None, gen_min_coverage=None, input_vcf=None,
                     process_input_positions_only=None, use_input_allele_only=None, use_sse_basecaller=None,
                     resolve_clipped_bases=None, prediction_precision=None, shift_likelihood_penalty=None,
                     minimum_sigma_prior=None, slope_sigma_prior=None, sigma_prior_weight=None, k_zero=None,
                     sse_relative_safety_level=None, tune_sbias=None, snp_min_coverage=None,
                     snp_min_cov_each_strand=None, snp_min_variant_score=None, snp_strand_bias=None,
                     snp_strand_bias_pval=None, snp_min_allele_freq=None, mnp_min_coverage=None,
                     mnp_min_cov_each_strand=None, mnp_min_variant_score=None, mnp_strand_bias=None,
                     mnp_strand_bias_pval=None, mnp_min_allele_freq=None, indel_min_coverage=None,
                     indel_min_cov_each_strand=None, indel_min_variant_score=None, indel_strand_bias=None,
                     indel_strand_bias_pval=None, indel_min_allele_freq=None, hotspot_min_coverage=None,
                     hotspot_min_cov_each_strand=None, hotspot_min_variant_score=None, hotspot_strand_bias=None,
                     hotspot_strand_bias_pval=None, hotspot_min_allele_freq=None, hp_max_length=None,
                     error_motifs=None, sse_prob_threshold=None, min_ratio_reads_non_sse_strand=None,
                     data_quality_stringency=None, read_rejection_threshold=None, filter_unusual_predictions=None,
                     filter_deletion_predictions=None, filter_insertion_predictions=None, heal_snps=None,
                     min_delta_for_flow=None, max_flows_to_test=None, outlier_probability=None, heavy_tailed=None,
                     suppress_recalibration=None, do_snp_realignment=None, do_mnp_realignment=None):

        options = ""
        options += " --input-bam %s" % input_bam
        options += " --reference %s" % reference_fasta
        options += " --output-dir %s" % out_dir
        options += " --output-vcf %s" % out_vcf

        options += " --sample-name %s" % sample_name if sample_name else ""
        options += " --num-variants-per-thread %s" % str(num_variants_per_thread) if num_variants_per_thread else ""
        options += " --parameters-file %s" % str(parameters_file) if parameters_file else ""
        options += " --force-sample-name %s" % str(force_sample_name) if force_sample_name else ""
        options += " --target-file %s" % str(target_file) if target_file else ""
        options += " --trim-ampliseq-primers %s" % str(trim_ampliseq_primers) if trim_ampliseq_primers else ""
        options += " --downsample-to-coverage %s" % str(downsample_to_coverage) if downsample_to_coverage else ""
        options += " --model-file %s" % str(model_file) if model_file else ""
        options += " --recal-model-hp-thres %s" % str(recal_model_hp_thres) if recal_model_hp_thres else ""
        options += " --output-vcf %s" % str(output_vcf) if output_vcf else ""
        options += " --suppress-reference-genotypes %s" % str(suppress_reference_genotypes) if suppress_reference_genotypes else ""
        options += " --suppress-no-calls %s" % str(suppress_no_calls) if suppress_no_calls else ""
        options += " --suppress-nocall-genotypes %s" % str(suppress_nocall_genotypes) if suppress_nocall_genotypes else ""
        options += " --allow-snps %s" % str(allow_snps) if allow_snps else ""
        options += " --allow-indels %s" % str(allow_indels) if allow_indels else ""
        options += " --allow-mnps %s" % str(allow_mnps) if allow_mnps else ""
        options += " --allow-complex %s" % str(allow_complex) if allow_complex else ""
        options += " --max-complex-gap %s" % str(max_complex_gap) if max_complex_gap else ""
        options += " --use-best-n-alleles %s" % str(use_best_n_alleles) if use_best_n_alleles else ""
        options += " --min-mapping-qv %s" % str(min_mapping_qv) if min_mapping_qv else ""
        options += " --read-snp-limit %s" % str(read_snp_limit) if read_snp_limit else ""
        options += " --read-max-mismatch-fraction %s" % str(read_max_mismatch_fraction) if read_max_mismatch_fraction else ""
        options += " --gen-min-alt-allele-freq %s" % str(gen_min_alt_allele_freq) if gen_min_alt_allele_freq else ""
        options += " --gen-min-indel-alt-allele-freq %s" % str(gen_min_indel_alt_allele_freq) if gen_min_indel_alt_allele_freq else ""
        options += " --gen-min-coverage %s" % str(gen_min_coverage) if gen_min_coverage else ""
        options += " --input-vcf %s" % str(input_vcf) if input_vcf else ""
        options += " --process-input-positions-only %s" % str(process_input_positions_only) if process_input_positions_only else ""
        options += " --use-input-allele-only %s" % str(use_input_allele_only) if use_input_allele_only else ""
        options += " --min-delta-for-flow %s" % str(min_delta_for_flow) if min_delta_for_flow else ""
        options += " --max-flows-to-test %s" % str(max_flows_to_test) if max_flows_to_test else ""
        options += " --outlier-probability %s" % str(outlier_probability) if outlier_probability else ""
        options += " --heavy-tailed %s" % str(heavy_tailed) if heavy_tailed else ""
        options += " --suppress-recalibration %s" % str(suppress_recalibration) if suppress_recalibration else ""
        options += " --do-snp-realignment %s" % str(do_snp_realignment) if do_snp_realignment else ""
        options += " --do-mnp-realignment %s" % str(do_mnp_realignment) if do_mnp_realignment else ""
        options += " --use-sse-basecaller %s" % str(use_sse_basecaller) if use_sse_basecaller else ""
        options += " --resolve-clipped-bases %s" % str(resolve_clipped_bases) if resolve_clipped_bases else ""
        options += " --prediction-precision %s" % str(prediction_precision) if prediction_precision else ""
        options += " --shift-likelihood-penalty %s" % str(shift_likelihood_penalty) if shift_likelihood_penalty else ""
        options += " --minimum-sigma-prior %s" % str(minimum_sigma_prior) if minimum_sigma_prior else ""
        options += " --slope-sigma-prior %s" % str(slope_sigma_prior) if slope_sigma_prior else ""
        options += " --sigma-prior-weight %s" % str(sigma_prior_weight) if sigma_prior_weight else ""
        options += " --k-zero %s" % str(k_zero) if k_zero else ""
        options += " --sse-relative-safety-level %s" % str(sse_relative_safety_level) if sse_relative_safety_level else ""
        options += " --tune-sbias %s" % str(tune_sbias) if tune_sbias else ""
        options += " --snp-min-coverage %s" % str(snp_min_coverage) if snp_min_coverage else ""
        options += " --snp-min-cov-each-strand %s" % str(snp_min_cov_each_strand) if snp_min_cov_each_strand else ""
        options += " --snp-min-variant-score %s" % str(snp_min_variant_score) if snp_min_variant_score else ""
        options += " --snp-strand-bias %s" % str(snp_strand_bias) if snp_strand_bias else ""
        options += " --snp-strand-bias-pval %s" % str(snp_strand_bias_pval) if snp_strand_bias_pval else ""
        options += " --snp-min-allele-freq %s" % str(snp_min_allele_freq) if snp_min_allele_freq else ""
        options += " --mnp-min-coverage %s" % str(mnp_min_coverage) if mnp_min_coverage else ""
        options += " --mnp-min-cov-each-strand %s" % str(mnp_min_cov_each_strand) if mnp_min_cov_each_strand else ""
        options += " --mnp-min-variant-score %s" % str(mnp_min_variant_score) if mnp_min_variant_score else ""
        options += " --mnp-strand-bias %s" % str(mnp_strand_bias) if mnp_strand_bias else ""
        options += " --mnp-strand-bias-pval %s" % str(mnp_strand_bias_pval) if mnp_strand_bias_pval else ""
        options += " --mnp-min-allele-freq %s" % str(mnp_min_allele_freq) if mnp_min_allele_freq else ""
        options += " --indel-min-coverage %s" % str(indel_min_coverage) if indel_min_coverage else ""
        options += " --indel-min-cov-each-strand %s" % str(indel_min_cov_each_strand) if indel_min_cov_each_strand else ""
        options += " --indel-min-variant-score %s" % str(indel_min_variant_score) if indel_min_variant_score else ""
        options += " --indel-strand-bias %s" % str(indel_strand_bias) if indel_strand_bias else ""
        options += " --indel-strand-bias-pval %s" % str(indel_strand_bias_pval) if indel_strand_bias_pval else ""
        options += " --indel-min-allele-freq %s" % str(indel_min_allele_freq) if indel_min_allele_freq else ""
        options += " --hotspot-min-coverage %s" % str(hotspot_min_coverage) if hotspot_min_coverage else ""
        options += " --hotspot-min-cov-each-strand %s" % str(hotspot_min_cov_each_strand) if hotspot_min_cov_each_strand else ""
        options += " --hotspot-min-variant-score %s" % str(hotspot_min_variant_score) if hotspot_min_variant_score else ""
        options += " --hotspot-strand-bias %s" % str(hotspot_strand_bias) if hotspot_strand_bias else ""
        options += " --hotspot-strand-bias-pval %s" % str(hotspot_strand_bias_pval) if hotspot_strand_bias_pval else ""
        options += " --hotspot-min-allele-freq %s" % str(hotspot_min_allele_freq) if hotspot_min_allele_freq else ""
        options += " --hp-max-length %s" % str(hp_max_length) if hp_max_length else ""
        options += " --error-motifs %s" % str(error_motifs) if error_motifs else ""
        options += " --sse-prob-threshold %s" % str(sse_prob_threshold) if sse_prob_threshold else ""
        options += " --min-ratio-reads-non-sse-strand %s" % str(min_ratio_reads_non_sse_strand) if min_ratio_reads_non_sse_strand else ""
        options += " --data-quality-stringency %s" % str(data_quality_stringency) if data_quality_stringency else ""
        options += " --read-rejection-threshold %s" % str(read_rejection_threshold) if read_rejection_threshold else ""
        options += " --filter-unusual-predictions %s" % str(filter_unusual_predictions) if filter_unusual_predictions else ""
        options += " --filter-deletion-predictions %s" % str(filter_deletion_predictions) if filter_deletion_predictions else ""
        options += " --filter-insertion-predictions %s" % str(filter_insertion_predictions) if filter_insertion_predictions else ""
        options += " --heal-snps %s" % str(heal_snps) if heal_snps else ""

        self.execute(options)

        """
        General options:

          -n,--num-threads                      INT         number of worker threads [2]
          -N,--num-variants-per-thread          INT         worker thread batch size [500]
             --parameters-file                  FILE        json file with algorithm control parameters [optional]

        Inputs:
          -r,--reference                        FILE        reference fasta file [required]
          -b,--input-bam                        FILE        bam file with mapped reads [required]
          -g,--sample-name                      STRING      sample for which variants are called (In case of input BAM files with multiple samples) [optional if there is only one sample]
             --force-sample-name                STRING      force all read groups to have this sample name [off]
          -t,--target-file                      FILE        only process targets in this bed file [optional]
             --trim-ampliseq-primers            on/off      match reads to targets and trim the ends that reach outside them [off]
          -D,--downsample-to-coverage           INT         ?? [2000]
             --model-file                       FILE        HP recalibration model input file.
             --recal-model-hp-thres             INT         Lower threshold for HP recalibration.

        Outputs:
          -O,--output-dir                       DIRECTORY   base directory for all output files [current dir]
          -o,--output-vcf                       FILE        vcf file with variant calling results [required]
             --suppress-reference-genotypes     on/off      write reference calls into the filtered variants vcf [on]
             --suppress-no-calls                on/off      write filtered variants into the filtered variants vcf [on]
             --suppress-nocall-genotypes        on/off      do not report a genotype for filtered variants [on]

        Variant candidate generation (FreeBayes):
             --allow-snps                       on/off      allow generation of snp candidates [on]
             --allow-indels                     on/off      allow generation of indel candidates [on]
             --allow-mnps                       on/off      allow generation of mnp candidates [on]
             --allow-complex                    on/off      allow generation of block substitution candidates [off]
             --max-complex-gap                  INT         maximum number of reference alleles between two alternate alleles to allow merging of the alternate alleles [1]
          -m,--use-best-n-alleles               INT         maximum number of snp alleles [2]
          -M,--min-mapping-qv                   INT         do not use reads with mapping quality below this [4]
          -U,--read-snp-limit                   INT         do not use reads with number of snps above this [10]
          -z,--read-max-mismatch-fraction       FLOAT       do not use reads with fraction of mismatches above this [1.0]
             --gen-min-alt-allele-freq          FLOAT       minimum required alt allele frequency to generate a candidate [0.2]
             --gen-min-indel-alt-allele-freq    FLOAT       minimum required alt allele frequency to generate a homopolymer indel candidate [0.2]
             --gen-min-coverage                 INT         minimum required coverage to generate a candidate [6]

        External variant candidates:
          -c,--input-vcf                        FILE        vcf.gz file (+.tbi) with additional candidate variant locations and alleles [optional]
             --process-input-positions-only     on/off      only generate candidates at locations from input-vcf [off]
             --use-input-allele-only            on/off      only consider provided alleles for locations in input-vcf [off]

        Variant candidate scoring options:
             --min-delta-for-flow               FLOAT       minimum prediction delta for scoring flows [0.1]
             --max-flows-to-test                INT         maximum number of scoring flows [10]
             --outlier-probability              FLOAT       probability for outlier reads [0.01]
             --heavy-tailed                     INT         degrees of freedom in t-dist modeling signal residual heavy tail [3]
             --suppress-recalibration           on/off      Suppress homopolymer recalibration [on].
             --do-snp-realignment               on/off      Realign reads in the vicinity of candidate snp variants [on].
             --do-mnp-realignment               on/off      Realign reads in the vicinity of candidate mnp variants [do-snp-realignment].

        Advanced variant candidate scoring options:
             --use-sse-basecaller               on/off      Switch to use the vectorized version of the basecaller [on].
             --resolve-clipped-bases            on/off      If 'true', the basecaller is used to solve soft clipped bases [off].
             --prediction-precision             FLOAT       prior weight in bias estimator [30.0]
             --shift-likelihood-penalty         FLOAT       penalize log-likelihood for solutions involving large systematic bias [0.3]
             --minimum-sigma-prior              FLOAT       prior variance per data point, constant [0.085]
             --slope-sigma-prior                FLOAT       prior rate of increase of variance over minimum by signal [0.0084]
             --sigma-prior-weight               FLOAT       weight of prior estimate of variance compared to observations [1.0]
             --k-zero                           FLOAT       variance increase for adding systematic bias [3.0]
             --sse-relative-safety-level        FLOAT       dampen strand bias detection for SSE events for low coverage [0.025]
             --tune-sbias                       FLOAT       dampen strand bias detection for low coverage [0.01]

        Variant filtering:
          -k,--snp-min-coverage                 INT         filter out snps with total coverage below this [6]
          -C,--snp-min-cov-each-strand          INT         filter out snps with coverage on either strand below this [1]
          -B,--snp-min-variant-score            FLOAT       filter out snps with QUAL score below this [2.5]
          -s,--snp-strand-bias                  FLOAT       filter out snps with strand bias above this [0.95] given pval < snp-strand-bias-pval
             --snp-strand-bias-pval             FLOAT       filter out snps with pval below this [1.0] given strand bias > snp-strand-bias
          -A,--snp-min-allele-freq              FLOAT       minimum required alt allele frequency for non-reference snp calls [0.2]
             --mnp-min-coverage                 INT         filter out mnps with total coverage below this [snp-min-coverage]
             --mnp-min-cov-each-strand          INT         filter out mnps with coverage on either strand below this, [snp-min-coverage-each-strand]
             --mnp-min-variant-score            FLOAT       filter out mnps with QUAL score below this [snp-min-variant-score]
             --mnp-strand-bias                  FLOAT       filter out mnps with strand bias above this [snp-strand-bias] given pval < mnp-strand-bias-pval
             --mnp-strand-bias-pval             FLOAT       filter out mnps with pval below this [snp-strand-bias-pval] given strand bias > mnp-strand-bias
             --mnp-min-allele-freq              FLOAT       minimum required alt allele frequency for non-reference mnp calls [snp-min-allele-freq]
             --indel-min-coverage               INT         filter out indels with total coverage below this [30]
             --indel-min-cov-each-strand        INT         filter out indels with coverage on either strand below this [1]
             --indel-min-variant-score          FLOAT       filter out indels with QUAL score below this [2.5]
          -S,--indel-strand-bias                FLOAT       filter out indels with strand bias above this [0.95] given pval < indel-strand-bias-pval
             --indel-strand-bias-pval           FLOAT       filter out indels with pval below this [1.0] given strand bias > indel-strand-bias
             --indel-min-allele-freq            FLOAT       minimum required alt allele frequency for non-reference indel call [0.2]
             --hotspot-min-coverage             INT         filter out hotspot variants with total coverage below this [6]
             --hotspot-min-cov-each-strand      INT         filter out hotspot variants with coverage on either strand below this [1]
             --hotspot-min-variant-score        FLOAT       filter out hotspot variants with QUAL score below this [2.5]
             --hotspot-strand-bias              FLOAT       filter out hotspot variants with strand bias above this [0.95] given pval < hotspot-strand-bias-pval
             --hotspot-strand-bias-pval         FLOAT       filter out hotspot variants with pval below this [1.0] given strand bias > hotspot-strand-bias
          -H,--hotspot-min-allele-freq          FLOAT       minimum required alt allele frequency for non-reference hotspot variant call [0.2]
          -L,--hp-max-length                    INT         filter out indels in homopolymers above this [8]
          -e,--error-motifs                     FILE        table of systematic error motifs and their error rates [optional]
             --sse-prob-threshold               FLOAT       filter out variants in motifs with error rates above this [0.2]
             --min-ratio-reads-non-sse-strand   FLOAT       minimum required alt allele frequency for variants with error motifs on opposite strand [0.2]
             --data-quality-stringency          FLOAT       minimum mean log-likelihood delta per read [4.0]
             --read-rejection-threshold         FLOAT       filter variants where large numbers of reads are rejected as outliers [0.5]
             --filter-unusual-predictions       FLOAT       posterior log likelihood threshold for accepting bias estimate [0.3]
             --filter-deletion-predictions      FLOAT       check post-evaluation systematic bias in deletions [100.0]
             --filter-insertion-predictions     FLOAT       check post-evaluation systematic bias in insertions [100.0]
             --heal-snps
        """


if __name__ == "__main__":

    input_bam = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/csa/alignment/001/Alsu24mc/tmap/uncorrected/001_trimmed_sorted_rm_pcr_chrom.bam"
    reference_fasta = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    out_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/csa/alignment/001/Alsu24mc/tmap/uncorrected/tvc"
    out_vcf = "out.vcf"
    TVC = TVC()
    #TVC.variant_call(input_bam, reference_fasta, out_vcf, out_dir=out_dir)

    TVC.variant_caller_pipeline(input_bam, reference_fasta, tvc_root_dir="/home/mahajrod/Repositories/genetic/NGS_tools/tmap/tvc-4.2.3-build", )
