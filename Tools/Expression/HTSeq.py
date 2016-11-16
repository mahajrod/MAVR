#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class HTSeq(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "htseq-count", path=path, max_threads=max_threads)

    def count(self, alignment_file, gff_file, output_file, samtype="bam", order=None, stranded_rnaseq="yes",
              min_alignment_quality=10, feature_type="exon", feature_id_attribute="gene_id", mode="union",
              suppress_progres_report=False):
        """
          -h, --help            show this help message and exit
          -f SAMTYPE, --format=SAMTYPE
                                type of <alignment_file> data, either 'sam' or 'bam'
                                (default: sam)
          -r ORDER, --order=ORDER
                                'pos' or 'name'. Sorting order of <alignment_file>
                                (default: name). Paired-end sequencing data must be
                                sorted either by position or by read name, and the
                                sorting order must be specified. Ignored for single-
                                end data.
          -s STRANDED, --stranded=STRANDED
                                whether the data is from a strand-specific assay.
                                Specify 'yes', 'no', or 'reverse' (default: yes).
                                'reverse' means 'yes' with reversed strand
                                interpretation
          -a MINAQUAL, --minaqual=MINAQUAL
                                skip all reads with alignment quality lower than the
                                given minimum value (default: 10)
          -t FEATURETYPE, --type=FEATURETYPE
                                feature type (3rd column in GFF file) to be used, all
                                features of other type are ignored (default, suitable
                                for Ensembl GTF files: exon)
          -i IDATTR, --idattr=IDATTR
                                GFF attribute to be used as feature ID (default,
                                suitable for Ensembl GTF files: gene_id)
          -m MODE, --mode=MODE  mode to handle reads overlapping more than one feature
                                (choices: union, intersection-strict, intersection-
                                nonempty; default: union)
          -o SAMOUT, --samout=SAMOUT
                                write out all SAM alignment records into an output SAM
                                file called SAMOUT, annotating each line with its
                                feature assignment (as an optional field with tag
                                'XF')
          -q, --quiet           suppress progress report
        """

        options = " -f %s" % samtype
        options += " -r %s" % order if order else ""
        options += " -s %s" % stranded_rnaseq if stranded_rnaseq else ""
        options += " -a %i" % min_alignment_quality if min_alignment_quality else ""
        options += " -t %s" % feature_type if feature_type else ""
        options += " -i %s" % feature_id_attribute if feature_id_attribute else ""
        options += " -m %s" % mode if mode else ""
        options += " -q" if suppress_progres_report else ""

        options += " %s" % alignment_file
        options += " %s" % gff_file
        options += " > %s" % output_file

        self.execute(options=options)
