__author__ = 'mahajrod'


from Tools.Abstract import Tool


class Tophat(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "tophat", path=path, max_threads=max_threads)

    def align(self, index, left_reads, right_reads=None, out_dir="tophat_output", read_mismatches=None,
              read_gap_length=None, read_edit_dist=None, read_realign_edit_dist=None, min_anchor=None,
              splice_mismatches=None, min_intron_length=None, max_intron_length=None, max_multihits=None,
              transcriptome_max_hits=None, max_insertion_length=None, max_deletion_length=None, quality=None,
              library_type=None, known_transcripts_file=None, transcriptome_index=None, raw_juncs_file=None,
              insertions_file=None, deletions_file=None, mate_inner_dist=None, mate_std_dev=None, no_novel_juncs=None,
              no_novel_indels=None, no_gtf_juncs=None, no_coverage_search=None, coverage_search=None,
              microexon_search=None, keep_tmp=None, tmp_dir=None, report_secondary_alignments=None, no_discordant=None,
              no_mixed=None, segment_mismatches=None, segment_length=None, bowtie_n=None, min_coverage_intron=None,
              max_coverage_intron=None, min_segment_intron=None, max_segment_intron=None, no_sort_bam=None,
              no_convert_bam=None, keep_fasta_order=None, allow_partial_mapping=None, bowtie2_mode=None,
              fusion_search=None, fusion_anchor_length=None, fusion_min_dist=None, fusion_read_mismatches=None,
              fusion_multireads=None, fusion_multipairs=None, fusion_ignore_chromosomes=None, rg_id=None,
              rg_sample=None, rg_library=None, rg_description=None, rg_platform_unit=None, rg_center=None,
              rg_date=None, rg_platform=None):
    
        """
        Usage:
        tophat [options] <bowtie_index> <reads1[,reads2,...]> [reads1[,reads2,...]] \
                                        [quals1,[quals2,...]] [quals1[,quals2,...]]

        Options:
            -v/--version
            -o/--output-dir                <string>    [ default: ./tophat_out         ]
            --bowtie1                                  [ default: bowtie2              ]
            -N/--read-mismatches           <int>       [ default: 2                    ]
            --read-gap-length              <int>       [ default: 2                    ]
            --read-edit-dist               <int>       [ default: 2                    ]
            --read-realign-edit-dist       <int>       [ default: "read-edit-dist" + 1 ]
            -a/--min-anchor                <int>       [ default: 8                    ]
            -m/--splice-mismatches         <0-2>       [ default: 0                    ]
            -i/--min-intron-length         <int>       [ default: 50                   ]
            -I/--max-intron-length         <int>       [ default: 500000               ]
            -g/--max-multihits             <int>       [ default: 20                   ]
            --suppress-hits
            -x/--transcriptome-max-hits    <int>       [ default: 60                   ]
            -M/--prefilter-multihits                   ( for -G/--GTF option, enable
                                                         an initial bowtie search
                                                         against the genome )
            --max-insertion-length         <int>       [ default: 3                    ]
            --max-deletion-length          <int>       [ default: 3                    ]
            --solexa-quals
            --solexa1.3-quals                          (same as phred64-quals)
            --phred64-quals                            (same as solexa1.3-quals)
            -Q/--quals
            --integer-quals
            -C/--color                                 (Solid - color space)
            --color-out
            --library-type                 <string>    (fr-unstranded, fr-firststrand,
                                                        fr-secondstrand)
            -p/--num-threads               <int>       [ default: 1                   ]
            -R/--resume                    <out_dir>   ( try to resume execution )
            -G/--GTF                       <filename>  (GTF/GFF with known transcripts)
            --transcriptome-index          <bwtidx>    (transcriptome bowtie index)
            -T/--transcriptome-only                    (map only to the transcriptome)
            -j/--raw-juncs                 <filename>
            --insertions                   <filename>
            --deletions                    <filename>
            -r/--mate-inner-dist           <int>       [ default: 50                  ]
            --mate-std-dev                 <int>       [ default: 20                  ]
            --no-novel-juncs
            --no-novel-indels
            --no-gtf-juncs
            --no-coverage-search
            --coverage-search
            --microexon-search
            --keep-tmp
            --tmp-dir                      <dirname>   [ default: <output_dir>/tmp ]
            -z/--zpacker                   <program>   [ default: gzip             ]
            -X/--unmapped-fifo                         [use mkfifo to compress more temporary
                                                         files for color space reads]

        Advanced Options:
            --report-secondary-alignments
            --no-discordant
            --no-mixed

            --segment-mismatches           <int>       [ default: 2                ]
            --segment-length               <int>       [ default: 25               ]

            --bowtie-n                                 [ default: bowtie -v        ]
            --min-coverage-intron          <int>       [ default: 50               ]
            --max-coverage-intron          <int>       [ default: 20000            ]
            --min-segment-intron           <int>       [ default: 50               ]
            --max-segment-intron           <int>       [ default: 500000           ]
            --no-sort-bam                              (Output BAM is not coordinate-sorted)
            --no-convert-bam                           (Do not output bam format.
                                                        Output is <output_dir>/accepted_hit.sam)
            --keep-fasta-order
            --allow-partial-mapping

        Bowtie2 related options:
          Preset options in --end-to-end mode (local alignment is not used in TopHat2)
            --b2-very-fast
            --b2-fast
            --b2-sensitive
            --b2-very-sensitive

          Alignment options
            --b2-N                         <int>       [ default: 0                ]
            --b2-L                         <int>       [ default: 20               ]
            --b2-i                         <func>      [ default: S,1,1.25         ]
            --b2-n-ceil                    <func>      [ default: L,0,0.15         ]
            --b2-gbar                      <int>       [ default: 4                ]

          Scoring options
            --b2-mp                        <int>,<int> [ default: 6,2              ]
            --b2-np                        <int>       [ default: 1                ]
            --b2-rdg                       <int>,<int> [ default: 5,3              ]
            --b2-rfg                       <int>,<int> [ default: 5,3              ]
            --b2-score-min                 <func>      [ default: L,-0.6,-0.6      ]

          Effort options
            --b2-D                         <int>       [ default: 15               ]
            --b2-R                         <int>       [ default: 2                ]

        Fusion related options:
            --fusion-search
            --fusion-anchor-length         <int>       [ default: 20               ]
            --fusion-min-dist              <int>       [ default: 10000000         ]
            --fusion-read-mismatches       <int>       [ default: 2                ]
            --fusion-multireads            <int>       [ default: 2                ]
            --fusion-multipairs            <int>       [ default: 2                ]
            --fusion-ignore-chromosomes    <list>      [ e.g, <chrM,chrX>          ]

            --fusion-do-not-resolve-conflicts          [this is for test purposes  ]

        SAM Header Options (for embedding sequencing run metadata in output):
            --rg-id                        <string>    (read group ID)
            --rg-sample                    <string>    (sample ID)
            --rg-library                   <string>    (library ID)
            --rg-description               <string>    (descriptive string, no tabs allowed)
            --rg-platform-unit             <string>    (e.g Illumina lane ID)
            --rg-center                    <string>    (sequencing center name)
            --rg-date                      <string>    (ISO 8601 date of the sequencing run)
            --rg-platform                  <string>    (Sequencing platform descriptor)

            for detailed help see http://tophat.cbcb.umd.edu/manual.html
        """
        options = ""
        options += " -o %s" % out_dir
        options += " --read-mismatches %i" % read_mismatches if read_mismatches else ""
        options += " --read-gap-length %i" % read_gap_length if read_gap_length else ""
        options += " --read-edit-dist %i" % read_edit_dist if read_edit_dist else ""
        options += " --read-realign-edit-dist %i" % read_realign_edit_dist if read_realign_edit_dist else ""
        options += " --min-anchor %i" % min_anchor if min_anchor else ""
        options += " --splice-mismatches %i" % splice_mismatches if splice_mismatches else ""
        options += " --min-intron-length %i" % min_intron_length if min_intron_length else ""
        options += " --max-intron-length %i" % max_intron_length if max_intron_length else ""
        options += " --max-multihits %i" % max_multihits if max_multihits else ""

        options += " --transcriptome-max-hits %i" % transcriptome_max_hits if transcriptome_max_hits else ""
        options += " --max-insertion-length %i" % max_insertion_length if max_insertion_length else ""
        options += " --max-deletion-length %i" % max_deletion_length if max_deletion_length else ""

        options += " --phred64" if quality == "phred64" else ""
        options += " --num-threads %i" % self.threads
        options += " --library-type %s" % library_type if library_type else ""
        options += " --GTF %s" % known_transcripts_file if known_transcripts_file else ""
        options += " --transcriptome-index %s" % transcriptome_index if transcriptome_index else ""
        options += " --raw-juncs %s" % raw_juncs_file if raw_juncs_file else ""
        options += " --insertions %s" % insertions_file if insertions_file else ""
        options += " --deletions %s" % deletions_file if deletions_file else ""
        options += " --mate-inner-dist %i" % mate_inner_dist if mate_inner_dist else ""
        options += " --mate-std-dev %i" % mate_std_dev if mate_std_dev else ""
        options += " --no-novel-juncs" if no_novel_juncs else ""
        options += " --no-novel-indels" if no_novel_indels else ""
        options += " --no-gtf-juncs" if no_gtf_juncs else ""
        options += " --no-coverage-search" if no_coverage_search else ""
        options += " --coverage-search" if coverage_search else ""
        options += " --microexon-search" if microexon_search else ""
        options += " --keep-tmp" if keep_tmp else ""
        options += " --tmp-dir %s" % tmp_dir if tmp_dir else ""

        options += " --report-secondary-alignments" if report_secondary_alignments else ""
        options += " --no-discordant" if no_discordant else ""
        options += " --no-mixed" if no_mixed else ""

        options += " --segment-mismatches %i" % segment_mismatches if segment_mismatches else ""
        options += " --segment-length %i" % segment_length if segment_length else ""
        options += " --bowtie-n" if bowtie_n else ""
        options += " --min-coverage-intron %i" % min_coverage_intron if min_coverage_intron else ""
        options += " --max-coverage-intron %i" % max_coverage_intron if max_coverage_intron else ""
        options += " --min-segment-intron %i" % min_segment_intron if min_segment_intron else ""
        options += " --max-segment-intron %i" % max_segment_intron if max_segment_intron else ""
        options += " --no-sort-bam" if no_sort_bam else ""
        options += " --no-convert-bam" if no_convert_bam else ""
        options += " --keep-fasta-order" if keep_fasta_order else ""
        options += " --allow-partial-mapping" if allow_partial_mapping else ""

        options += " --b2-%s" % bowtie2_mode if bowtie2_mode else ""
        #allowed mode: very-fast, fast, sensitive, very-sensitive

        options += " --fusion-search" if fusion_search else ""
        options += " --fusion-anchor-length %i" % fusion_anchor_length if fusion_anchor_length else ""
        options += " --fusion-min-dist %i" % fusion_min_dist if fusion_min_dist else ""
        options += " --fusion-read-mismatches %i" % fusion_read_mismatches if fusion_read_mismatches else ""
        options += " --fusion-multireads %i" % fusion_multireads if fusion_multireads else ""
        options += " --fusion-multipairs %i" % fusion_multipairs if fusion_multipairs else ""
        options += " --fusion-ignore-chromosomes %s" % fusion_ignore_chromosomes if fusion_ignore_chromosomes else ""

        options += " --rg-id %s" % rg_id if rg_id else ""
        options += " --rg-sample %s" % rg_sample if rg_sample else ""
        options += " --rg-library %s" % rg_library if rg_library else ""
        options += " --rg-description %s" % rg_description if rg_description else ""
        options += " --rg-platform-unit %s" % rg_platform_unit if rg_platform_unit else ""
        options += " --rg-center %s" % rg_center if rg_center else ""
        options += " --rg-date %s" % rg_date if rg_date else ""
        options += " --rg-platform %s" % rg_platform if rg_platform else ""

        options += " %s" % index
        options += " %s" % left_reads
        options += " %s" % right_reads if right_reads else ""

        self.execute(options, cmd="tophat")