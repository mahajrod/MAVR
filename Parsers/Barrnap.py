#!/usr/bin/env python

from collections import OrderedDict

from Bio import  SeqRecord, SeqFeature

from CustomCollections.GeneralCollections import TwoLvlDict


class RecordBARRNAP():
    def __init__(self, start, end, strand, type, length, expected_length, partial):

        self.start = start                                          # int
        self.end = end                                              # int
        self.strand = strand                                        # str
        self.type = type                                            # str
        self.length = length                                        # int
        self.expected_length = expected_length                      # int
        self.partial = partial                                      # boolean
        self.length_ratio = float(length) / float(expected_length) if expected_length else None

    def attributes_string(self):
        attributes_string = "Len=%i;Exp_len=%i;Partial=%s" % (self.length, self.expected_length, str(self.partial))

        attributes_string += (";Ratio=%f" % self.length_ratio) if self.length_ratio else ""

        return attributes_string

    def gff_str(self):
        #source	type	start	end	score	strand	phase	attributes
        return "Barrnap\t%s\t%i\t%i\t.\t%s\t.\t%s" % (self.type, self.start, self.end, self.strand,
                                                      self.attributes_string())

"""
    def simple_gff_str(self):
        # seqid	source	type	start	end	score	strand	phase	attributes
        attributes_string = "Period=%i;N_copies=%.1f;Pattern=%s" % (self.period, self.number_of_copies, self.pattern)

        return "TRF\trepeat\t%i\t%i\t.\t.\t.\t%s" % (self.start, self.end, attributes_string)

    def table_str_short(self):

        return "%i\t%i\t%i\t%.1f\t%s\t%s" % (self.start, self.end, self.period, self.number_of_copies,
                                             self.pattern, self.tandem_repeat)

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

"""


class CollectionBARRNAP():

    def __init__(self, record_dict=None, e_value_cutoff=None, partial_threshold=None, reject_threshold=None,
                 barrnap_file=None, from_file=True):
        if from_file:
            self.records = OrderedDict()
            with open(barrnap_file, "r") as fd:
                chrom = None
                for line in fd:
                    if "Setting evalue cutoff" in line:
                        self.e_value_cutoff = line.strip().split()[-1]
                    elif "Will tag genes" in line:
                        self.partial_threshold = float(line.strip().split()[-4])
                    elif "Will reject genes" in line:
                        self.reject_threshold = float(line.strip().split()[-4])
                    if "Found:" in line:
                        self._add_record(line)
        else:
            self.records = record_dict
            self.e_value_cutoff = e_value_cutoff            # e-value cuttof for hit
            self.partial_threshold = partial_threshold      # length threshold to mark gene as partial (if below)
            self.reject_threshold = reject_threshold        # length threshold to reject gene (if below)

    def _add_record(self, line):
        line_list = line.strip().split("Found: ")[-1].split()
        type = line_list[0]
        chrom = line_list[1]
        length, expected_length = map(int, line_list[2].split("=")[-1].split("/"))
        start, end = map(int, line_list[3].split(".."))
        strand = line_list[4]
        partial = "partial" in line_list[-1]
        record = RecordBARRNAP(start, end, strand, type, length, expected_length, partial)
        if chrom not in self.records:
            self.records[chrom] = []
        self.records[chrom].append(record)

    def get_annotated_types(self):

        annotated_types = set()
        for chrom in self.records:
            for record in self.records[chrom]:
                annotated_types.add(record.type)

        return annotated_types

    def count_types(self, output_file=None, total_output_file=None, return_mode="chrom"):

        annotated_types = self.get_annotated_types()
        count_dict = TwoLvlDict()
        total_count_dict = OrderedDict()

        for type in annotated_types:
            total_count_dict[type] = OrderedDict()
            total_count_dict[type]["complete"] = 0
            total_count_dict[type]["partial"] = 0

        for chrom in self.records:
            count_dict[chrom] = OrderedDict()
            for type in annotated_types:
                count_dict[chrom][type] = 0

        for chrom in self.records:
            for record in self.records[chrom]:
                count_dict[chrom][record.type] += 1
                if record.partial:
                    total_count_dict[record.type]["partial"] += 1
                else:
                    total_count_dict[record.type]["complete"] += 1

        if output_file:
            count_dict.write(output_file)

        if total_output_file:
            with open(total_output_file, "w") as out_fd:
                out_fd.write("#rRNA\tComplete\tPartial%s\n" % ("(<%.2f of expected length)" if
                                                               self.partial_threshold else ""))
                for type in total_count_dict:
                    out_fd.write("%s\t%i\t%i\n" % (type, total_count_dict[type]["complete"],
                                                   total_count_dict[type]["partial"]))

        if return_mode == "chrom":
            return count_dict
        elif return_mode == "total":
            return total_count_dict
        elif return_mode == "both":
            return count_dict, total_count_dict
        else:
            raise ValueError("Unknown return type. Allowed variants: 'chrom', 'total', 'both'")

    """
    def count_types(self, output_file=None):

        annotated_types = self.get_annotated_types()
        total_count_dict = OrderedDict()

        for type in annotated_types:
            total_count_dict[type] = OrderedDict()
            total_count_dict[type]["complete"] = 0
            total_count_dict[type]["partial"] = 0

        for chrom in self.records:
            for record in self.records[chrom]:
                if record.partial:
                    total_count_dict[record.type]["partial"] += 1
                else:
                    total_count_dict[record.type]["complete"] += 1
        if output_file:
            with open(output_file) as out_fd:
                out_fd.write("#rRNA\tNumber\n")
                for type in total_count_dict:
                    out_fd.write("%s\t%i\n" % (type, total_count_dict[type]))

        return total_count_dict
    """

    def write_gff(self, out_file, write_chr_string=False):
        with open(out_file, "w") as out_fd:
            for chrom in self.records:
                if write_chr_string:
                    out_fd.write("#%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.gff_str()))


    def write(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#>chrom\n#\tstart\tend\tperiod\tnumber_of_copies\tpattern\ttandem_repeat\n")
            for chrom in self.records:
                out_fd.write(">%s\n" % chrom)
                for record in self.records[chrom]:
                    out_fd.write("\t%s\n" % record.table_str_short())

    def write_short_table(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#chrom\tstart\tend\tperiod\tnumber_of_copies\tpattern\ttandem_repeat\n")
            for chrom in self.records:
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.table_str_short()))

    def write_wide_table(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write("#chrom\tstart\tend\tperiod\tnumber_of_copies\tconsensus_pattern_size\tpercent_of_matches\tpercent_of_indels\talignment_score\tnucleotide_composition\tentropy\tpattern\ttandem_repeat")
            for chrom in self.records:
                for record in self.records[chrom]:
                    out_fd.write("%s\t%s\n" % (chrom, record.table_str_long()))

