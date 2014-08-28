
from Tools.Abstract import Tool


class BWA(Tool):
    def __init__(self, path=""):
        super(Tool).__init__("bwa", path=path)

    def index(self, fasta_file, prefix=None):
        index_prefix = ""
        if prefix:
            index_prefix = "-p %s" % prefix
        options = "%s %s" % (index_prefix, fasta_file)
        self.execute(options, command="bwa index")
        #os.system("bwa index %s %s" % (index_prefix, fasta_file))

    def mem(self, index, forward_reads, reverse_reads=None, output_file="alignment.sam", max_threads=4):
        reads = forward_reads
        if reverse_reads:
            reads += " %s" % reverse_reads
        options = "-t %i %s %s > %s" % (max_threads, index, reads, output_file)
        self.execute(options, command="bwa mem")
        #os.system("bwa mem -t %i %s %s > %s" % (max_threads, index, reads, output_file))