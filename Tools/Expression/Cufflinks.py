#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class Cufflinks(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cufflinks", path=path, max_threads=max_threads)

    def assembly_transcripts(self, alignment, output_dir=None, GTF=None, GTF_guide=None, mask_file=None,
                             frag_bias_correct=None, label=None, library_type=None, library_norm_method=None, seed=None,
                             frag_len_mean=None, frag_len_std_dev=None, max_mle_iterations=None,
                             num_frag_count_draws=None, num_frag_assign_draws=None, max_frag_multihits=None,
                             max_intron_length=None, min_frags_per_transfrag=None, overhang_tolerance=None,
                             max_bundle_length=None, max_bundle_frags=None, min_intron_length=None,
                             trim_3_avgcov_thresh=None, overlap_radius=None, three_overhang_tolerance=None,
                             intron_overhang_tolerance=None, min_isoform_fraction=None, pre_mrna_fraction=None,
                             junc_alpha=None, small_anchor_fraction=None, trim_3_dropoff_frac=None,
                             max_multiread_fraction=None, multi_read_correct=None, no_effective_length_correction=None,
                             no_length_correction=None, compatible_hits_norm=None, total_hits_norm=None,
                             no_faux_reads=None, verbose=None, quiet=None, no_update_check=None):

        options = " --num-threads %i" % self.threads
        options += " --output-dir %s" % output_dir if output_dir else ""
        options += " --GTF %s" % GTF if GTF else ""
        options += " --GTF-guide %s" % GTF_guide if GTF_guide else ""
        options += " --mask-file %s" % mask_file if mask_file else ""
        options += " --frag-bias-correct %s" % frag_bias_correct if frag_bias_correct else ""
        options += " --label %s" % label if label else ""
        options += " --library-type %s" % library_type if library_type else ""
        options += " --library-norm-method %s" % library_norm_method if library_norm_method else ""

        options += " --seed %i" % seed if seed else ""
        options += " --frag-len-mean %i" % frag_len_mean if frag_len_mean else ""
        options += " --frag-len-std-dev %i" % frag_len_std_dev if frag_len_std_dev else ""
        options += " --max-mle-iterations %i" % max_mle_iterations if max_mle_iterations else ""

        options += " --num-frag-count-draws %i" % num_frag_count_draws if num_frag_count_draws else ""
        options += " --num-frag-assign-draws %i" % num_frag_assign_draws if num_frag_assign_draws else ""
        options += " --max-frag-multihits %i" % max_frag_multihits if max_frag_multihits else ""
        options += " --max-intron-length %i" % max_intron_length if max_intron_length else ""
        options += " --min-frags-per-transfrag %i" % min_frags_per_transfrag if min_frags_per_transfrag else ""
        options += " --overhang-tolerance %i" % overhang_tolerance if overhang_tolerance else ""
        options += " --max-bundle-length %i" % max_bundle_length if max_bundle_length else ""
        options += " --max-bundle-frags %i" % max_bundle_frags if max_bundle_frags else ""
        options += " --min-intron-length %i" % min_intron_length if min_intron_length else ""
        options += " --trim-3-avgcov-thresh %i" % trim_3_avgcov_thresh if trim_3_avgcov_thresh else ""
        options += " --overlap-radius %i" % overlap_radius if overlap_radius else ""

        options += " --3-overhang-tolerance %i" % three_overhang_tolerance if three_overhang_tolerance else ""
        options += " --intron-overhang-tolerance %i" % intron_overhang_tolerance if intron_overhang_tolerance else ""

        options += " --min-isoform-fraction %f" % min_isoform_fraction if min_isoform_fraction else ""
        options += " --pre-mrna-fraction %f" % pre_mrna_fraction if pre_mrna_fraction else ""
        options += " --junc-alpha %f" % junc_alpha if junc_alpha else ""
        options += " --small-anchor-fraction %f" % small_anchor_fraction if small_anchor_fraction else ""
        options += " --trim-3-dropoff-frac %f" % trim_3_dropoff_frac if trim_3_dropoff_frac else ""
        options += " --max-multiread-fraction %f" % max_multiread_fraction if max_multiread_fraction else ""

        options += " --multi-read-correct" if multi_read_correct else ""
        options += " --no-effective-length-correction" if no_effective_length_correction else ""
        options += " --no-length-correction" if no_length_correction else ""
        options += " --compatible-hits-norm" if compatible_hits_norm else ""
        options += " --total-hits-norm" if total_hits_norm else ""
        options += " --no-faux-reads" if no_faux_reads else ""
        options += " --verbose" if verbose else ""
        options += " --quiet" if quiet else ""
        options += " --no-update-check" if no_update_check else ""

        options += " %s" % alignment

        self.execute(options=options)


class Cuffmerge(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cuffmerge", path=path, max_threads=max_threads)

    def merge(self, gtf_list_file, output_dir="./", reference_gtf=None, reference_fasta=None, min_isoform_fraction=None, keep_tmp=None):

        options = " --num-threads %i" % self.threads
        options += " -o %s" % output_dir
        options += " --ref-gtf %s " % reference_gtf if reference_gtf else ""
        options += " --ref_sequence %s" % reference_fasta if reference_fasta else ""
        options += " --min-isoform-fraction %f" % min_isoform_fraction if min_isoform_fraction else ""
        options += " --keep-tmp" if keep_tmp else ""
        options += " %s" % gtf_list_file

        self.execute(options=options)


class Cuffquant(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cuffquant", path=path, max_threads=max_threads)

    def quant(self, output_dir, gtf_file, sam_files, mask_file=None, frag_bias_correct=None, library_type=None, seed=None, quiet=None,
              verbose=None, frag_len_mean=None, frag_len_std_dev=None, max_mle_iterations=None, no_update_check=None,
              max_bundle_frags=None, max_frag_multihits=None, no_effective_length_correction=None,
              no_length_correction=None, no_scv_correction=None, read_skip_fraction=None, no_read_pairs=None,
              min_alignment_count=None):

        options = " --num-threads %i" % self.threads
        options += " --output-dir %s" % output_dir if output_dir else ""
        options += " --mask-file %s" % mask_file if mask_file else ""
        options += " --frag-bias-correct %s" % frag_bias_correct if frag_bias_correct else ""
        options += " --library-type %s" % library_type if library_type else ""
        options += " --seed %i" % seed if seed else ""
        options += " --frag-len-mean %i" % frag_len_mean if frag_len_mean else ""
        options += " --frag-len-std-dev %i" % frag_len_std_dev if frag_len_std_dev else ""
        options += " --max-mle-iterations %i" % max_mle_iterations if max_mle_iterations else ""
        options += " --max-frag-multihits %i" % max_frag_multihits if max_frag_multihits else ""
        options += " --max-bundle-frags %i" % max_bundle_frags if max_bundle_frags else ""
        options += " --min-alignment-count %i" % min_alignment_count if min_alignment_count else ""
        options += " --read-skip-fraction %f" % read_skip_fraction if read_skip_fraction else ""
        options += " --no-effective-length-correction" if no_effective_length_correction else ""
        options += " --no-length-correction" if no_length_correction else ""
        options += " --no-read-pairs" if no_read_pairs else ""
        options += " --no-scv-correction" if no_scv_correction else ""
        options += " --verbose" if verbose else ""
        options += " --quiet" if quiet else ""
        options += " --no-update-check" if no_update_check else ""

        options += " %s" % gtf_file
        options += " %s" % (sam_files if isinstance(sam_files, str) else " ".join(sam_files))

        self.execute(options=options)


class Cuffnorm(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cuffnorm", path=path, max_threads=max_threads)

    def normalize(self, output_dir, gtf_file, cxb_files_two_lvl_list, labels_list=None, norm_standards_file=None,
                  library_norm_method=None, library_type=None, seed=None, quiet=None,
                  verbose=None, no_update_check=None, output_format=None,
                  compatible_hits_norm=None, total_hits_norm=None):

        options = " --num-threads %i" % self.threads
        options += " --output-dir %s" % output_dir if output_dir else ""
        options += " --output_format" % output_format if output_format else ""
        options += " --library-type %s" % library_type if library_type else ""
        options += " --norm_standards_file %s" % norm_standards_file if norm_standards_file else ""
        options += " --library-norm-method %s" % library_norm_method if library_norm_method else ""
        options += (" --labels %s" % ",".join(labels_list)) if labels_list else ""
        options += " --seed %i" % seed if seed else ""
        options += " --compatible-hits-norm" if compatible_hits_norm else ""
        options += " --total-hits-norm" if total_hits_norm else ""
        options += " --verbose" if verbose else ""
        options += " --no-update-check" if no_update_check else ""
        options += " --quiet" if quiet else ""
        options += " %s" % gtf_file
        options += " %s" % " ".join(map(lambda s: s if isinstance(s, str) else ",".join(s), cxb_files_two_lvl_list))

        self.execute(options=options)


class Cuffdiff(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cuffdiff", path=path, max_threads=max_threads)

    def compare(self, gtf_file, sam_files_two_lvl_list, output_dir=None, mask_file=None, frag_bias_correct=None,
                library_type=None, library_norm_method=None, seed=None, frag_len_mean=None, frag_len_std_dev=None,
                max_mle_iterations=None, num_frag_count_draws=None, num_frag_assign_draws=None, max_frag_multihits=None,
                max_bundle_frags=None, multi_read_correct=None, no_effective_length_correction=None,
                no_length_correction=None, compatible_hits_norm=None, total_hits_norm=None, verbose=None, quiet=None,
                no_update_check=None, labels_list=None, dispersion_method=None, min_reps_for_js_test=None,
                no_scv_correction=None, min_alignment_count=None, read_skip_fraction=None, no_read_pairs=None,
                FDR=None, no_diff=None, no_js_tests=None, time_series=None
                ):

        options = " --num-threads %i" % self.threads
        options += " --output-dir %s" % output_dir if output_dir else ""
        options += (" --labels %s" % ",".join(labels_list)) if labels_list else ""
        options += " --mask-file %s" % mask_file if mask_file else ""
        options += " --multi-read-correct" if multi_read_correct else ""
        options += " --library-type %s" % library_type if library_type else ""
        options += " --dispersion-method %s" % dispersion_method if dispersion_method else ""

        options += " --library-norm-method %s" % library_norm_method if library_norm_method else ""
        options += " --frag-len-mean %i" % frag_len_mean if frag_len_mean else ""
        options += " --frag-len-std-dev %i" % frag_len_std_dev if frag_len_std_dev else ""

        options += " --max-bundle-frags %i" % max_bundle_frags if max_bundle_frags else ""
        options += " --min-reps-for-js-test %i" % min_reps_for_js_test if min_reps_for_js_test else ""

        options += " --compatible-hits-norm" if compatible_hits_norm else ""
        options += " --total-hits-norm" if total_hits_norm else ""

        options += " --no-scv-correction" if no_scv_correction else ""
        options += " --verbose" if verbose else ""
        options += " --quiet" if quiet else ""

        options += " --min-alignment-count %i" % min_alignment_count if min_alignment_count else ""
        options += " --read-skip-fraction %f" % read_skip_fraction if read_skip_fraction else ""
        options += " --no-effective-length-correction" if no_effective_length_correction else ""
        options += " --no-length-correction" if no_length_correction else ""

        options += " --num-frag-count-draws %i" % num_frag_count_draws if num_frag_count_draws else ""
        options += " --num-frag-assign-draws %i" % num_frag_assign_draws if num_frag_assign_draws else ""
        options += " --max-frag-multihits %i" % max_frag_multihits if max_frag_multihits else ""
        options += " --no-update-check" if no_update_check else ""
        options += " --seed %i" % seed if seed else ""
        options += " --max-mle-iterations %i" % max_mle_iterations if max_mle_iterations else ""
        options += " --no-read-pairs" if no_read_pairs else ""
        options += " --frag-bias-correct %s" % frag_bias_correct if frag_bias_correct else ""

        options += " --FDR %f" % FDR if FDR else ""
        options += " --no-diff" if no_diff else ""
        options += " --no-js-tests" if no_js_tests else ""
        options += " --time-series" if time_series else ""

        options += " %s" % gtf_file
        options += " %s" % " ".join(map(lambda s: s if isinstance(s, str) else ",".join(s), sam_files_two_lvl_list))

        self.execute(options=options)


class Cuffcompare(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cuffcompare", path=path, max_threads=max_threads)

    def compare(self):

        options = ""
        self.execute(options=options)
        pass


class Gtf2sam(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "gtf_to_sam", path=path, max_threads=max_threads)

    def convert(self, gtf_list, output_sam, reference_fasta=None, use_fpkm_instead_isoform_fraction=False):

        options = " --reference-seq %s" % reference_fasta if reference_fasta else ""
        options += " --raw-fpkm" if use_fpkm_instead_isoform_fraction else ""
        options += " %s" % (gtf_list if isinstance(gtf_list, str) else ",".join(gtf_list))
        options += " %s" % output_sam


class Gffread(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "gffread", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(input_gff_file, genomic_fasta_file=None, output_cds_file=None, output_protein_file=None,
                             output_transcripts_file=None, coding_only=False):

        options = " %s" % input_gff_file
        options += " -C" if coding_only else ""
        options += " -g %s" % genomic_fasta_file if genomic_fasta_file else ""
        options += " -x %s" % output_cds_file if output_cds_file else ""
        options += " -w %s" % output_transcripts_file if output_transcripts_file else ""
        options += " -y %s" % output_protein_file if output_protein_file else ""
        return options

    def extract_transcript_sequences(self, input_gff_file, genomic_fasta_file, output_prefix, coding_only=False):

        output_cds_file = "%s.cds" % output_prefix
        output_transcripts_file = "%s.transcript" % output_prefix
        output_protein_file = "%s.protein" % output_prefix
        options = self.parse_common_options(input_gff_file,
                                            coding_only=coding_only,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_cds_file=output_cds_file,
                                            output_protein_file=output_protein_file,
                                            output_transcripts_file=output_transcripts_file)
        #print options
        self.execute(options=options)

    def extract_cds(self, input_gff_file, genomic_fasta_file, output_file):

        options = self.parse_common_options(input_gff_file,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_cds_file=output_file,
                                            coding_only=True,
                                            output_transcripts_file=None)
        #print options
        self.execute(options=options)

    def extract_mRNA(self, input_gff_file, genomic_fasta_file, output_file, coding_only=False):

        options = self.parse_common_options(input_gff_file,
                                            coding_only=coding_only,
                                            genomic_fasta_file=genomic_fasta_file,
                                            output_transcripts_file=output_file)
        #print options
        self.execute(options=options)

if __name__ == "__main__":
    pass
