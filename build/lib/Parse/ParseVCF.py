#!/usr/bin/python2


class MutRecord(object):
    """docstring for MutRecord"""
    def __init__(self, chrom=None, pos=None, identificator=None, ref=None, alt=None, qual=None,
                 filter_list=[], info_dict={}, format_list=[], sample_list=[], strand=None):
        super(MutRecord, self).__init__()
        self.chrom = chrom
        self.pos = pos
        self.id = identificator
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_list
        self.info = info_dict
        self.format = format_list
        self.samples = sample_list
        self.strand = strand
        self.frequency = None
        self.distance = None

    def calc_freq(self):
        self.frequency = 0
        for sample in self.samples:
            gentype = sample.data.GT.split("|")
            self.frequency += int(gentype[0]) + int(gentype[1])


def parse_vcf(vcf_filename):
    fd = open(vcf_filename, "w")
    format = fd.readline().strip()[2:].split("=")
    source = fd.readline().strip()[2:].split("=")
    reference = fd.readline().strip()[2:].split("=")

    rest_of_metadata = []
    for line in fd:
        if line[:6] == "#CHROM":
            columns_names = line.strip()[1:].split("\t")
            break
        rest_of_metadata.append(line.strip())
    if "FORMAT" in columns_names:
        samples_names = columns_names[9:]
        info_list = columns_names[7].split(";")
        for line in fd:
            splited = line.strip().split("\t")
            yield MutRecord(chrom=splited[0],
                            pos=splited[1],
                            identificator=splited[2],
                            ref=splited[3],
                            alt=splited[4],
                            qual=splited[5],
                            filter_list=splited[6],
                            )

    fd.close()