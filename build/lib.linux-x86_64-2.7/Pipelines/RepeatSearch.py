import os
from Tools.BLAST import make_blast_plus_db


def sliding_window_blastn(sequence, blast_db, window_size=100, step=50, blast_task="blastn-short", threads=8, output_file="output.bl"):
    seq_length = len(sequence)
    start = 0
    with open(output_file, "w"):
        pass
    blast_string = "blastn -task %s -num_threads %i -db %s -num_descriptions 500000000 -num_alignments 500000000  -word_size 4 -evalue 1e300 -ungapped -outfmt \"10 sacc qstart qend qlen sstart send sstrand length evalue\"  >> %s" \
                   % (blast_task, threads, blast_db, output_file)

    while (seq_length - start) >= 2 * window_size:
        os.system("echo \"%s\" |" % (str(sequence[start:start+window_size])) + blast_string)
        start += step

    os.system("echo \"%s\" |" % (str(sequence[start])) + blast_string)


if __name__ == "__main__":
    pass