
from Tools.Abstract import Tool


class Cufflinks(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cufflinks", path=path, max_threads=max_threads)
        
    def assembly(self, alignment, out_dir="cufflinks_output", seed=None, GTF=None, GTF_guide=None, mask_file=None,
                 frag_bias_correct=None, multi_read_correct=None, library_type=None, library_norm_method=None,
                 frag_len_mean=None, frag_len_std_dev=None, max_mle_iterations=None, compatible_hits_norm=None,
                 total_hits_norm=None, num_frag_count_draws=None, num_frag_assign_draws=None, max_frag_multihits=None,
                 no_effective_length_correction=None, no_length_correction=None, upper_quartile_norm=None,
                 label=None, min_isoform_fraction=None, pre_mrna_fraction=None, max_intron_length=None,
                 junc_alpha=None, small_anchor_fraction=None, min_frags_per_transfrag=None, overhang_tolerance=None,
                 max_bundle_length=None, max_bundle_frags=None, min_intron_length=None, trim_3_avgcov_thresh=None,
                 trim_3_dropoff_frac=None, max_multiread_fraction=None, overlap_radius=None, no_faux_reads=None,
                 overhang_tolerance_3=None, intron_overhang_tolerance=None):
        """
        Usage:   cufflinks [options] <hits.sam>
        General Options:
          -o/--output-dir              write all output files to this directory              [ default:     ./ ]
          -p/--num-threads             number of threads used during analysis                [ default:      1 ]
          --seed                       value of random number generator seed                 [ default:      0 ]
          -G/--GTF                     quantitate against reference transcript annotations                      
          -g/--GTF-guide               use reference transcript annotation to guide assembly                   
          -M/--mask-file               ignore all alignment within transcripts in this file                     
          -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
          -u/--multi-read-correct      use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]
          --library-type               library prep used for input reads                     [ default:  below ]
          --library-norm-method        Method used to normalize library sizes                [ default:  below ]
        
        Advanced Abundance Estimation Options:
          -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]
          -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]
          --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]
          --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:  FALSE ]
          --total-hits-norm            count all hits for normalization                      [ default:  TRUE  ]
          --num-frag-count-draws       Number of fragment generation samples                 [ default:    100 ]
          --num-frag-assign-draws      Number of fragment assignment samples per generation  [ default:     50 ]
          --max-frag-multihits         Maximum number of alignments allowed per fragment     [ default: unlim  ]
          --no-effective-length-correction   No effective length correction                  [ default:  FALSE ]
          --no-length-correction       No length correction                                  [ default:  FALSE ]
          -N/--upper-quartile-norm     Deprecated, use --library-norm-method                 [    DEPRECATED   ]
          --raw-mapped-norm            Deprecated, use --library-norm-method                 [    DEPRECATED   ]
        
        Advanced Assembly Options:
          -L/--label                   assembled transcripts have this ID prefix             [ default:   CUFF ]
          -F/--min-isoform-fraction    suppress transcripts below this abundance level       [ default:   0.10 ]
          -j/--pre-mrna-fraction       suppress intra-intronic transcripts below this level  [ default:   0.15 ]
          -I/--max-intron-length       ignore alignments with gaps longer than this          [ default: 300000 ]
          -a/--junc-alpha              alpha for junction binomial test filter               [ default:  0.001 ]
          -A/--small-anchor-fraction   percent read overhang taken as 'suspiciously small'   [ default:   0.09 ]
          --min-frags-per-transfrag    minimum number of fragments needed for new transfrags [ default:     10 ]
          --overhang-tolerance         number of terminal exon bp to tolerate in introns     [ default:      8 ]
          --max-bundle-length          maximum genomic length allowed for a given bundle     [ default:3500000 ]
          --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]
          --min-intron-length          minimum intron size allowed in genome                 [ default:     50 ]
          --trim-3-avgcov-thresh       minimum avg coverage required to attempt 3' trimming  [ default:     10 ]
          --trim-3-dropoff-frac        fraction of avg coverage below which to trim 3' end   [ default:    0.1 ]
          --max-multiread-fraction     maximum fraction of allowed multireads per transcript [ default:   0.75 ]
          --overlap-radius             maximum gap size to fill between transfrags (in bp)   [ default:     50 ]
        
        Advanced Reference Annotation Guided Assembly Options:
          --no-faux-reads              disable tiling by faux reads                          [ default:  FALSE ]
          --3-overhang-tolerance       overhang allowed on 3' end when merging with reference[ default:    600 ]
          --intron-overhang-tolerance  overhang allowed inside reference intron when merging [ default:     30 ]
        
        Advanced Program Behavior Options:
          -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
          -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
          --no-update-check            do not contact server to check for update availability[ default:  FALSE ]
        
        Supported library types:
            ff-firststrand
            ff-secondstrand
            ff-unstranded
            fr-firststrand
            fr-secondstrand
            fr-unstranded (default)
            transfrags
        
        Supported library normalization methods:
            classic-fpkm
            poisson
        """
        options = ""
        options += " --output-dir %s" % out_dir
        options += " --num-threads %i" % self.threads
        options += " --seed %i" % seed if seed else ""
        options += " --GTF %s" % GTF if GTF else ""
        options += " --GTF-guide %s" % GTF_guide if GTF_guide else ""
        options += " --mask-file %s" % mask_file if mask_file else ""
        options += " --frag-bias-correct %s" % frag_bias_correct if frag_bias_correct else ""
        options += " --multi-read-correct" if multi_read_correct else ""
        options += " --library-type %s" % library_type if library_type else ""
        # allowed library types:
        # ff-firststrand
        # ff-secondstrand
        # ff-unstranded
        # fr-firststrand
        # fr-secondstrand
        # fr-unstranded (default)
        # transfrags

        options += " --library-norm-method %s" % library_norm_method if library_norm_method else ""
        # allowed library normalizationmethods:
        # classic-fpkm
	    # poisson

        options += " --frag-len-mean %i" % frag_len_mean if frag_len_mean else ""
        options += " --frag-len-std-dev %i" % frag_len_std_dev if frag_len_std_dev else ""
        options += " --max-mle-iterations %i" % max_mle_iterations if max_mle_iterations else ""
        options += " --compatible-hits-norm" if compatible_hits_norm else ""
        options += " --total-hits-norm" if total_hits_norm else ""
        options += " --num-frag-count-draws %i" % num_frag_count_draws if num_frag_count_draws else ""
        options += " --num-frag-assign-draws %i" % num_frag_assign_draws if num_frag_assign_draws else ""
        options += " --max-frag-multihits %i" % max_frag_multihits if max_frag_multihits else ""
        options += " --no-effective-length-correction" if no_effective_length_correction else ""
        options += " --no-length-correction" if no_length_correction else ""
        options += " --upper-quartile-norm" if upper_quartile_norm else ""

        options += " --label %s" % label if label else ""
        options += " --min-isoform-fraction %f" % min_isoform_fraction if min_isoform_fraction else ""
        options += " --pre-mrna-fraction %f" % pre_mrna_fraction if pre_mrna_fraction else ""
        options += " --max-intron-length %i" % max_intron_length if max_intron_length else ""
        options += " --junc-alpha %f" % junc_alpha if junc_alpha else ""
        options += " --small-anchor-fraction %f" % small_anchor_fraction if small_anchor_fraction else ""
        options += " --min-frags-per-transfrag %i" % min_frags_per_transfrag if min_frags_per_transfrag else ""
        options += " --overhang-tolerance %i" % overhang_tolerance if overhang_tolerance else ""
        options += " --max-bundle-length %i" % max_bundle_length if max_bundle_length else ""
        options += " --max-bundle-frags %i" % max_bundle_frags if max_bundle_frags else ""
        options += " --min-intron-length %i" % min_intron_length if min_intron_length else ""
        options += " --trim-3-avgcov-thresh %i" % trim_3_avgcov_thresh if trim_3_avgcov_thresh else ""
        options += " --trim-3-dropoff-frac %f" % trim_3_dropoff_frac if trim_3_dropoff_frac else ""
        options += " --max-multiread-fraction %f" % max_multiread_fraction if max_multiread_fraction else ""
        options += " --overlap-radius %i" % overlap_radius if overlap_radius else ""
        options += " --no-faux-reads" if no_faux_reads else ""
        options += " --3-overhang-tolerance %i" % overhang_tolerance_3 if overhang_tolerance_3 else ""
        options += " --intron-overhang-tolerance %i" % intron_overhang_tolerance if intron_overhang_tolerance else ""

        options += " %s" % alignment

        self.execute(options, cmd="cufflinks")