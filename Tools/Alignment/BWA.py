
from Tools.Abstract import Tool


class BWA(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "bwa", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

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

    def align(self, index,
              forward_reads_list=None,
              reverse_reads_list=None,
              unpaired_reads_list=None,
              quality_score="phred33",
              output_prefix="alignment",
              output_format="bam",
              read_group_name="reads",
              PU="x",
              SM="sample",
              platform="Illumina",
              LB="x",
              sort_by_coordinate=True,
              sort_by_name=False,
              max_per_sorting_thread_memory=None
              ):
        """
        Common method for all aligners
        """

        if (not forward_reads_list) and (not reverse_reads_list):
            raise ValueError("No read files were given")

        reads = forward_reads_list if isinstance(forward_reads_list, str) else  forward_reads_list[0]
        if reverse_reads_list:
            reads += " %s" % reverse_reads_list if isinstance(reverse_reads_list, str) else reverse_reads_list[0]

        options = " -t %i" % self.threads
        options += " %s" % index
        options += " %s " % reads
        options += " -R \'@RG\\tID:%s\\tPU:%s\\tSM:%s\\tPL:%s\\tLB:%s\'" % (read_group_name, PU, SM, platform, LB)

        if sort_by_coordinate or sort_by_name:
            if sort_by_coordinate and sort_by_name:
                raise ValueError("Sorting by both coordinate and read name was requested")
            options += " | samtools sort"
            if sort_by_name:
                options += " -n"
            options += " -@ %i" % self.threads
            options += " -m %s" % max_per_sorting_thread_memory if max_per_sorting_thread_memory else self.max_per_thread_memory
            options += " -O %s" % output_format.upper()

        options += " > %s.%s" % (output_prefix, output_format)

        self.execute(options, cmd="bwa mem")


    def aln(self, in_db_fasta, in_query_fastq, out_sai, max_edit_dist=0.04, max_gaps=1, max_gap_ext=-1, disallow_lond_del_limit=16,
            disallow_indel_limit=5, max_seed_edit_dist=2, mismatch_penalty=3, gap_open_penalty=11,
            gap_extension_penalty=4, trim=0, barcode_length=0):

        options = "-t %i -n %f -o %i -e %i -d %i -i %s -k %i -M %i -O %i -E %i -q %i -B %i" %\
                  (self.threads, max_edit_dist, max_gaps, max_gap_ext, disallow_lond_del_limit, disallow_indel_limit,
                   max_seed_edit_dist, mismatch_penalty, gap_open_penalty, gap_extension_penalty, trim, barcode_length)

        options += " %s %s %s" % (in_db_fasta, in_query_fastq, out_sai)

        self.execute(options, cmd="bwa aln")