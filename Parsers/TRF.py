#!/usr/bin/env python


class RecordTRF():
    def __init__(self, chrom, start, end, period, number_of_copies, consensus_pattern_size,
                 percent_of_matches, percent_of_indels, alignment_score, nucleotide_percent_list, entropy,
                 pattern=None, tandem_repeat=None):
        # TRF repeat string description:
        # Indices of the repeat relative to the start of the sequence.
        # Period size of the repeat.
        # Number of copies aligned with the consensus pattern.
        # Size of consensus pattern (may differ slightly from the period size).
        # Percent of matches between adjacent copies overall.
        # Percent of indels between adjacent copies overall.
        # Alignment score.
        # Percent composition for each of the four nucleotides.
        # Entropy measure based on percent composition.
        self.chrom = chrom                                          # str
        self.start = start                                          # int
        self.end = end                                              # int
        self.period = period                                        # int
        self.number_of_copies = number_of_copies                    # float
        self.consensus_pattern_size = consensus_pattern_size        # int
        self.percent_of_matches = percent_of_matches                # int
        self.percent_of_indels = percent_of_indels                  # int
        self.alignment_score = alignment_score                      # int
        self.nucleotide_percent_list = nucleotide_percent_list      # list of four ints
        self.entropy = entropy                                      # float
        self.pattern = pattern                                      # str or None
        self.tandem_repeat = tandem_repeat                          # str or None

    def gff_str(self):
        # seqid	source	type	start	end	score	strand	phase	attributes
        attributes_string = "Period=%i;N_copies=%.1f;Cons_pat_size=%i;Pers_matches=%i;Pers_indels=%i;Align_score=%i" \
                            % (self.period, self.number_of_copies, self.consensus_pattern_size,
                               self.percent_of_matches, self.percent_of_indels, self.alignment_score)
        nuc_composition = ",".join(map(lambda x: str(x), self.nucleotide_percent_list))

        attributes_string += ";Nuc_composition=%s;Entropy=%.2f" % (nuc_composition, self.entropy)
        return "%s\tTRF\trepeat\t%i\t%i\t.\t.\t.\t%s" % (self.chrom, self.start, self.end, attributes_string)


class CollectionTRF():

    def __init__(self, parameters=None, record_list=None, trf_file=None, from_file=True):
        self.linkage_dict = None
        if from_file:
            self.records = []
            with open(trf_file, "r") as fd:
                tmp = next(fd)
                while tmp[0:8] != "Sequence":
                        tmp = next(fd)
                while True:
                    #line = fd.readline()
                    while tmp[0:8] != "Sequence":
                        tmp = next(fd)
                        #print tmp
                    chrom = tmp.strip().split()[1]
                    while tmp[0:10] != "Parameters":
                        tmp = next(fd)
                        #print tmp
                    print tmp
                    self.parameters = list(map(lambda x: int(x), tmp.strip().split()[1:]))
                    tmp = next(fd)
                    while tmp == "\n":
                        tmp = next(fd)
                    print tmp
                    while tmp != "\n" and tmp != "":
                        self._add_record(tmp, chrom)
                        tmp = next(fd)
        else:
            self.records = record_list
            self.parameters = parameters

    def _add_record(self, line, chrom):
        line_list = line.strip().split()
        start, end, period = list(map(lambda x: int(x), line_list[0:3]))
        consensus_pattern_size, percent_of_matches, percent_of_indels, alignment_score = \
            list(map(lambda x: int(x), line_list[4:8]))
        nucleotide_percent_list = list(map(lambda x: int(x), line_list[8:-3]))
        record = RecordTRF(chrom, start, end, period, float(line_list[3]), consensus_pattern_size,
                           percent_of_matches, percent_of_indels, alignment_score,
                           nucleotide_percent_list,
                           float(line_list[-3]),
                           pattern=line_list[-2], tandem_repeat=line_list[-1])
        self.records.append(record)

    def write_gff(self, out_file):
        with open(out_file, "w") as out_fd:
            for record in self.records:
                out_fd.write(record.gff_str() + "\n")


if __name__ == "__main__":
    trf_coll = CollectionTRF(trf_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masking/TRF/LAN210_v0.10m_masked_repeatmasker.fasta.2.7.7.80.10.50.500.dat", from_file=True)
    trf_coll.write_gff("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masking/TRF/trf.gff")