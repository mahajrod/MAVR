
from Tools.Abstract import Tool


class BWA(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bwa", path=path, max_threads=max_threads)

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

        options = "-t %i -n %f -o %i -e %i -d %i -i %s -k %i -M %i -O %i -E %i -q %i -B %i" %\
                  (self.threads, max_edit_dist, max_gaps, max_gap_ext, disallow_lond_del_limit, disallow_indel_limit,
                   max_seed_edit_dist, mismatch_penalty, gap_open_penalty, gap_extension_penalty, trim, barcode_length)

        options += " %s %s %s" % (in_db_fasta, in_query_fastq, out_sai)

        self.execute(options, cmd="bwa aln")