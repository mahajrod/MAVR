
from Tools.Abstract import Tool


class BWA(Tool):
    def __init__(self, path="", max_threads=4):
        super(Tool).__init__("bwa", path=path, max_threads=max_threads)

    def index(self, fasta_file, prefix=None):
        index_prefix = ""
        if prefix:
            index_prefix = "-p %s" % prefix
        options = "%s %s" % (index_prefix, fasta_file)
        self.execute(options, cmd="bwa index")
        #os.system("bwa index %s %s" % (index_prefix, fasta_file))

    def mem(self, index, forward_reads, reverse_reads=None, output_file="alignment.sam"):
        reads = forward_reads
        if reverse_reads:
            reads += " %s" % reverse_reads
        options = "-t %i %s %s > %s" % (self.threads, index, reads, output_file)
        self.execute(options, cmd="bwa mem")
        #os.system("bwa mem -t %i %s %s > %s" % (max_threads, index, reads, output_file))

    def aln(self, in_db_fasta, in_query_fastq, out_sai, max_edit_dist=0.04, max_gaps=1, max_gap_ext=-1, disallow_lond_del_limit=16,
            disallow_indel_limit=5, max_seed_edit_dist=2, mismatch_penalty=3, gap_open_penalty=11,
            gap_extension_penalty=4, trim=0, barcode_length=0):
        """
        -n NUM 	Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
        -o INT 	Maximum number of gap opens [1]
        -e INT 	Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
        -d INT 	Disallow a long deletion within INT bp towards the 3’-end [16]
        -i INT 	Disallow an indel within INT bp towards the ends [5]
        -l INT 	Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’. [inf]
        -k INT 	Maximum edit distance in the seed [2]
        -t INT 	Number of threads (multi-threading mode) [1]
        -M INT 	Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc). [3]
        -O INT 	Gap open penalty [11]
        -E INT 	Gap extension penalty [4]
        -B INT 	Length of barcode starting from the 5’-end. When INT is positive, the barcode of each read will be trimmed before mapping and will be written at the BC SAM tag. For paired-end reads, the barcode from both ends are concatenated. [0]
        -q INT 	Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]



        -R INT 	Proceed with suboptimal alignments if there are no more than INT equally best hits. This option only affects paired-end mapping. Increasing this threshold helps to improve the pairing accuracy at the cost of speed, especially for short reads (~32bp).
        -c 	Reverse query but not complement it, which is required for alignment in the color space. (Disabled since 0.6.x)
        -N 	Disable iterative search. All hits with no more than maxDiff differences will be found. This mode is much slower than the default.
        -I 	The input is in the Illumina 1.3+ read format (quality equals ASCII-64).
        -b 	Specify the input read sequence file is the BAM format. For paired-end data, two ends in a pair must be grouped together and options -1 or -2 are usually applied to specify which end should be mapped. Typical command lines for mapping pair-end data in the BAM format are:

        bwa aln ref.fa -b1 reads.bam > 1.sai
        bwa aln ref.fa -b2 reads.bam > 2.sai
        bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > aln.sam
        -0 	When -b is specified, only use single-end reads in mapping.
        -1 	When -b is specified, only use the first read in a read pair in mapping (skip single-end reads and the second reads).
        -2 	When -b is specified, only use the second read in a read pair in mapping.
        """
        options = "-t %i -n %f -o %i -e %i -d %i -i %s -k %i -M %i -O %i -E %i -q %i -B %i" %\
                  (self.threads, max_edit_dist, max_gaps, max_gap_ext, disallow_lond_del_limit, disallow_indel_limit,
                   max_seed_edit_dist, mismatch_penalty, gap_open_penalty, gap_extension_penalty, trim, barcode_length)

        options += " %s %s %s" % (in_db_fasta, in_query_fastq, out_sai)

        self.execute(options, cmd="bwa aln")