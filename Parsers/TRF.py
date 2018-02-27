#!/usr/bin/env python

from collections import OrderedDict


class RecordTRF():
    def __init__(self, start, end, period, number_of_copies, consensus_pattern_size,
                 percent_of_matches, percent_of_indels, alignment_score, nucleotide_percent_list, entropy,
                 pattern=None, tandem_repeat=None, record_id=None):
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
        self.id = record_id
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

    def attributes_string(self):
        attributes_string = "ID=%s;" % str(self.id) if self.id else ""
        attributes_string += "Period=%i;N_copies=%.1f;Pattern=%s;Cons_pat_size=%i;Pers_matches=%i;Pers_indels=%i;Align_score=%i" \
                             % (self.period, self.number_of_copies, self.pattern, self.consensus_pattern_size,
                                self.percent_of_matches, self.percent_of_indels, self.alignment_score)
        nuc_composition = ",".join(map(lambda x: str(x), self.nucleotide_percent_list))

        attributes_string += ";Nuc_composition=%s;Entropy=%.2f" % (nuc_composition, self.entropy)

        return attributes_string

    def atributes_string_with_rep_seq(self):
        return self.attributes_string() + ";seq=%s" % self.tandem_repeat

    def gff_str_with_rep_seq(self):
        return "TRF\trepeat\t%i\t%i\t.\t.\t.\t%s" % (self.start, self.end, self.atributes_string_with_rep_seq())

    def gff_str(self):
        # seqid	source	type	start	end	score	strand	phase	attributes
        return "TRF\trepeat\t%i\t%i\t.\t.\t.\t%s" % (self.start, self.end, self.attributes_string())

    def simple_gff_str(self):
        # seqid	source	type	start	end	score	strand	phase	attributes
        attributes_string = "ID=%s;" % self.id if self.id else ""
        attributes_string += "Period=%i;N_copies=%.1f;Pattern=%s" % (self.period, self.number_of_copies, self.pattern)

        return "TRF\trepeat\t%i\t%i\t.\t.\t.\t%s" % (self.start, self.end, attributes_string)

    def table_str_short(self):

        return "%i\t%i\t%i\t%.1f\t%s\t%s" % (self.start, self.end, self.period, self.number_of_copies,
                                             self.pattern, self.tandem_repeat)

    def simple_bed_str(self):

        return "%i\t%i\t%s\t.\t." % (self.start, self.end, self.id)

    def table_str_long(self):

        return "%i\t%i\t%i\t%.1f\t%i\t%i\t%i\t%i\t%s\t%.2f\t%s\t%s" % (self.start, self.end, self.period,
                                                                       self.number_of_copies,
                                                                       self.consensus_pattern_size,
                                                                       self.percent_of_matches,
                                                                       self.percent_of_indels,
                                                                       self.alignment_score,
                                                                       ",".join(map(lambda x: str(x), self.nucleotide_percent_list)),
                                                                       self.entropy,
                                                                       self.pattern, self.tandem_repeat)


class CollectionTRF():

    def __init__(self, parameters=None, record_list=None, trf_file=None, from_file=True, add_ids=True,
                 record_id_prefix="TRF"):
        self.linkage_dict = None
        record_index = 0
        if from_file:
            self.records = {}
            for trf_filename in [trf_file] if isinstance(trf_file, str) else trf_file:
                with open(trf_filename, "r") as fd:
                    chrom = None
                    for line in fd:
                        #tmp = line.strip()
                        if line[0:8] == "Sequence":
                            chrom = line.strip().split()[1]
                        elif not chrom:
                            continue
                        elif line[0:10] == "Parameters":
                            self.parameters = list(map(lambda x: int(x), line.strip().split()[1:]))
                        elif line != "\n":
                            self._add_record(line, chrom, record_id="%s%i" % (record_id_prefix, record_index) if add_ids else None)
                            record_index += 1
                """
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
                    while (tmp != "\n") and (tmp != "") and (tmp[0:8] != "Sequence"):
                        self._add_record(tmp, chrom)
                        try:
                         tmp = next(fd)
                        except StopIteration:
                            break
                """
        else:
            self.records = record_list
            self.parameters = parameters

    def _add_record(self, line, chrom, record_id=None):
        line_list = line.strip().split()
        start, end, period = list(map(lambda x: int(x), line_list[0:3]))
        consensus_pattern_size, percent_of_matches, percent_of_indels, alignment_score = \
            list(map(lambda x: int(x), line_list[4:8]))
        nucleotide_percent_list = list(map(lambda x: int(x), line_list[8:-3]))
        record = RecordTRF(start, end, period, float(line_list[3]), consensus_pattern_size,
                           percent_of_matches, percent_of_indels, alignment_score,
                           nucleotide_percent_list,
                           float(line_list[-3]),
                           pattern=line_list[-2], tandem_repeat=line_list[-1],
                           record_id=record_id)
        if chrom not in self.records:
            self.records[chrom] = []
        self.records[chrom].append(record)

    def write_gff_with_rep_seqs(self, out_file, write_chr_string=False):
        with open(out_file, "w") as out_fd:
            for chrom in self.records:
                if write_chr_string:
                    out_fd.write("#%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.gff_str_with_rep_seq()))

    def write_gff(self, out_file, write_chr_string=False):
        with open(out_file, "w") as out_fd:
            for chrom in self.records:
                if write_chr_string:
                    out_fd.write("#%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.gff_str()))

    def write_simple_gff(self, out_file, write_chr_string=False):
        with open(out_file, "w") as out_fd:
            for chrom in self.records:
                if write_chr_string:
                    out_fd.write("#%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.simple_gff_str()))

    def write(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#>chrom\n#\tstart\tend\tperiod\tnumber_of_copies\tpattern\ttandem_repeat\n")
            for chrom in self.records:
                out_fd.write(">%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("\t%s\n" % record.table_str_short())

    def write_bed(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#chrom\tstart\tend\tid\tscore\tstrand\n")
            for chrom in self.records:
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.simple_bed_str()))

    def write_short_table(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#chrom\tstart\tend\tperiod\tnumber_of_copies\tpattern\ttandem_repeat\n")
            for chrom in self.records:
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.table_str_short()))

    def write_wide_table(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#chrom\tstart\tend\tperiod\tnumber_of_copies\tconsensus_pattern_size\tpercent_of_matches\tpercent_of_indels\talignment_score\tnucleotide_composition\tentropy\tpattern\ttandem_repeat\n")
            for chrom in self.records:
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.table_str_long()))

    def write_fasta(self, out_file):
        with open(out_file, "w") as out_fd:
            for chrom in self.records:
                for record in self.records[chrom]:
                    if record.id and record.tandem_repeat:
                        out_fd.write(">%s location:%s,%i,%i pattern:%s period:%i copies:%f matches,%%:%i "
                                     "indels,%%:%i align_score:%i entropy:%f\n%s\n" % (record.id,
                                                                                       chrom,
                                                                                       record.start,
                                                                                       record.end,
                                                                                       record.pattern,
                                                                                       record.period,
                                                                                       record.number_of_copies,
                                                                                       record.percent_of_matches,
                                                                                       record.percent_of_indels,
                                                                                       record.alignment_score,
                                                                                       record.entropy,
                                                                                       record.tandem_repeat))

    def get_id_based_dict(self):
        """
        returns following dictionary {trf_id: [chrom, start, end]}
        """
        id_based_dict = OrderedDict()
        for chrom in self.records:
            for record in self.records[chrom]:
                id_based_dict[record.id] = [chrom, record.start, record.end]

        return id_based_dict


if __name__ == "__main__":
    trf_coll = CollectionTRF(trf_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masking/TRF/LAN210_v0.10m_masked_repeatmasker.fasta.2.7.7.80.10.50.500.dat", from_file=True)
    trf_coll.write_gff("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masking/TRF/trf.gff")
