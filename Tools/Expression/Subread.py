#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
from math import ceil
from CustomCollections.GeneralCollections import SynDict, TwoLvlDict
from Tools.Abstract import Tool


class Subread(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "featureCounts", path=path, max_threads=max_threads)

    def count(self, alignment_file, gff_file, output_file, count_multimapped_reads=False, annotation_file_type="GTF",
              min_read_fraction_overlap=None, feature_type_to_use=None, attribute_type_to_use=None, stranded=None):
        """
          Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...

## Mandatory arguments:

  -a <string>         Name of an annotation file. GTF/GFF format by default.
                      See -F option for more format information. Inbuilt
                      annotations (SAF format) is available in 'annotation'
                      directory of the package.

  -o <string>         Name of the output file including read counts. A separate
                      file including summary statistics of counting results is
                      also included in the output ('<string>.summary')

  input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be
                      either name or location sorted. If not files provided,
                      <stdin> input is expected.

## Optional arguments:
# Annotation

  -F <string>         Specify format of the provided annotation file. Acceptable
                      formats include 'GTF' (or compatible GFF format) and
                      'SAF'. 'GTF' by default.  For SAF format, please refer to
                      Users Guide.

  -t <string>         Specify feature type in GTF annotation. 'exon' by
                      default. Features used for read counting will be
                      extracted from annotation using the provided value.

  -g <string>         Specify attribute type in GTF annotation. 'gene_id' by
                      default. Meta-features used for read counting will be
                      extracted from annotation using the provided value.

  -A <string>         Provide a chromosome name alias file to match chr names in
                      annotation with those in the reads. This should be a two-
                      column comma-delimited text file. Its first column should
                      include chr names in the annotation and its second column
                      should include chr names in the reads. Chr names are case
                      sensitive. No column header should be included in the
                      file.

# Level of summarization

  -f                  Perform read counting at feature level (eg. counting
                      reads for exons rather than genes).

# Overlap between reads and features

  -O                  Assign reads to all their overlapping meta-features (or
                      features if -f is specified).

  --minOverlap <int>  Minimum number of overlapping bases in a read that is
                      required for read assignment. 1 by default. Number of
                      overlapping bases is counted from both reads if paired
                      end. If a negative value is provided, then a gap of up
                      to specified size will be allowed between read and the
                      feature that the read is assigned to.

  --fracOverlap <float> Minimum fraction of overlapping bases in a read that is
                      required for read assignment. Value should be within range
                      [0,1]. 0 by default. Number of overlapping bases is
                      counted from both reads if paired end. Both this option
                      and '--minOverlap' option need to be satisfied for read
                      assignment.

  --largestOverlap    Assign reads to a meta-feature/feature that has the
                      largest number of overlapping bases.

  --readExtension5 <int> Reads are extended upstream by <int> bases from their
                      5' end.

  --readExtension3 <int> Reads are extended upstream by <int> bases from their
                      3' end.

  --read2pos <5:3>    Reduce reads to their 5' most base or 3' most base. Read
                      counting is then performed based on the single base the
                      read is reduced to.

# Multi-mapping reads

  -M                  Multi-mapping reads will also be counted. For a multi-
                      mapping read, all its reported alignments will be
                      counted. The 'NH' tag in BAM/SAM input is used to detect
                      multi-mapping reads.

# Fractional counting

  --fraction          Assign fractional counts to features. This option must
                      be used together with '-M' or '-O' or both. When '-M' is
                      specified, each reported alignment from a multi-mapping
                      read (identified via 'NH' tag) will carry a fractional
                      count of 1/x, instead of 1 (one), where x is the total
                      number of alignments reported for the same read. When '-O'
                      is specified, each overlapping feature will receive a
                      fractional count of 1/y, where y is the total number of
                      features overlapping with the read. When both '-M' and
                      '-O' are specified, each alignment will carry a fractional
                      count of 1/(x*y).

# Read filtering

  -Q <int>            The minimum mapping quality score a read must satisfy in
                      order to be counted. For paired-end reads, at least one
                      end should satisfy this criteria. 0 by default.

  --splitOnly         Count split alignments only (ie. alignments with CIGAR
                      string containing 'N'). An example of split alignments is
                      exon-spanning reads in RNA-seq data.

  --nonSplitOnly      If specified, only non-split alignments (CIGAR strings do
                      not contain letter 'N') will be counted. All the other
                      alignments will be ignored.

  --primary           Count primary alignments only. Primary alignments are
                      identified using bit 0x100 in SAM/BAM FLAG field.

  --ignoreDup         Ignore duplicate reads in read counting. Duplicate reads
                      are identified using bit Ox400 in BAM/SAM FLAG field. The
                      whole read pair is ignored if one of the reads is a
                      duplicate read for paired end data.

# Strandness

  -s <int>            Perform strand-specific read counting. Acceptable values:
                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).
                      0 by default.

# Exon-exon junctions

  -J                  Count number of reads supporting each exon-exon junction.
                      Junctions were identified from those exon-spanning reads
                      in the input (containing 'N' in CIGAR string). Counting
                      results are saved to a file named '<output_file>.jcounts'

  -G <string>         Provide the name of a FASTA-format file that contains the
                      reference sequences used in read mapping that produced the
                      provided SAM/BAM files. This optional argument can be used
                      with '-J' option to improve read counting for junctions.

# Parameters specific to paired end reads

  -p                  If specified, fragments (or templates) will be counted
                      instead of reads. This option is only applicable for
                      paired-end reads.

  -B                  Only count read pairs that have both ends aligned.

  -P                  Check validity of paired-end distance when counting read
                      pairs. Use -d and -D to set thresholds.

  -d <int>            Minimum fragment/template length, 50 by default.

  -D <int>            Maximum fragment/template length, 600 by default.

  -C                  Do not count read pairs that have their two ends mapping
                      to different chromosomes or mapping to same chromosome
                      but on different strands.

  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from
                      the same pair are required to be located next to each
                      other in the input.

# Number of CPU threads

  -T <int>            Number of the threads. 1 by default.

# Miscellaneous

  -R                  Output detailed assignment result for each read. A text
                      file will be generated for each input file, including
                      names of reads and meta-features/features reads were
                      assigned to. See Users Guide for more details.

  --tmpDir <string>   Directory under which intermediate files are saved (later
                      removed). By default, intermediate files will be saved to
                      the directory specified in '-o' argument.

  --maxMOp <int>      Maximum number of 'M' operations allowed in a CIGAR
                      string. 10 by default. Both 'X' and '=' are treated as 'M'
                      and adjacent 'M' operations are merged in the CIGAR
                      string.

  -v                  Output version of the program.

        """

        options = " -T %i" % self.threads

        options += " -M" if count_multimapped_reads else ""
        options += " --fracOverlap %f" % min_read_fraction_overlap if min_read_fraction_overlap else ""
        options += " -F %s" % annotation_file_type

        options += " -t %s" % feature_type_to_use if feature_type_to_use else ""
        options += " -g %s" % attribute_type_to_use if attribute_type_to_use else ""
        options += " -F %s" % annotation_file_type

        options += " -a %s" % gff_file
        options += " -o %s" % output_file
        options += " %s" % alignment_file
        options += " -s %i" % stranded if stranded else ""

        self.execute(options=options)

    def count_miRNA_reads(self, alignment_file, gff_file, output_prefix, annotation_file_type="GTF",
                          min_read_fraction_overlap=1.0, feature_type_to_use=None, attribute_type_to_use=None,
                          sample_name=None, stranded=1):

        no_multimapped_read_counts = "%s.no_multimapped_reads.count" % output_prefix
        with_multimapped_read_counts = "%s.with_multimapped_reads.count" % output_prefix
        all_adjusted_read_counts = "%s.all_adjusted_reads.count" % output_prefix

        self.count(alignment_file, gff_file, no_multimapped_read_counts, annotation_file_type=annotation_file_type,
                   min_read_fraction_overlap=min_read_fraction_overlap, feature_type_to_use=feature_type_to_use,
                   attribute_type_to_use=attribute_type_to_use, stranded=stranded)

        self.count(alignment_file, gff_file, with_multimapped_read_counts, count_multimapped_reads=True,
                   annotation_file_type=annotation_file_type,
                   min_read_fraction_overlap=min_read_fraction_overlap, feature_type_to_use=feature_type_to_use,
                   attribute_type_to_use=attribute_type_to_use, stranded=stranded)

        no_multimapped_read_count_dict = SynDict(filename=no_multimapped_read_counts, comments_prefix="#",
                                                 key_index=0, value_index=6, expression=int, header=True)
        with_multimapped_read_count_dict = SynDict(filename=with_multimapped_read_counts, comments_prefix="#",
                                                   key_index=0, value_index=6, expression=int, header=True)
        similar_feature_number_dict = SynDict(filename=with_multimapped_read_counts, comments_prefix="#", header=True,
                                              key_index=0, value_index=1, expression=lambda s: len(s.split(";")))

        sample_nameeee = sample_name if sample_name else similar_feature_number_dict.header.split()[6]

        all_adjusted_read_count_dict = SynDict()
        all_adjusted_read_count_dict.header = ".\t%s" % sample_nameeee

        #print no_multimapped_read_count_dict
        #print with_multimapped_read_count_dict
        #print similar_feature_number_dict

        for feature_id in no_multimapped_read_count_dict:
            all_adjusted_read_count_dict[feature_id] = ceil(float(no_multimapped_read_count_dict[feature_id]) + \
                                                            (float(with_multimapped_read_count_dict[feature_id]) - float(no_multimapped_read_count_dict[feature_id])) / float(similar_feature_number_dict[feature_id]))

        all_adjusted_read_count_dict.write(all_adjusted_read_counts, header=True)






