#!/usr/bin/env python

from RouToolPa.Parsers.Abstract import Record, Metadata, Collection


class MetadataGFF(Metadata):
    # TODO: rewrite for general case
    pass


class RecordGFF(Record):

    def __init__(self, chrom, source, type, start, end, score, strand, quality, attributes):
        self.chrom = chrom                    # str
        self.source = source                    # str
        self.type = type                        # str
        self.start = start                      # int
        self.end = end                          # int
        self.score = score                      # float or '.'
        self.strand = strand                    # str : '+' or '-' or '.'
        self.quality = quality                  # float or '.'
        self.attributes = attributes            # generally dict, but maybe a str if unparsed

    def __str__(self):
        if isinstance(self.attributes, dict):
            attribute_string = ";".join(["%s=%s" % (key, self.attributes[key]) for key in self.attributes])
        else:
            attribute_string = self.attributes
        return "%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s" % (self.chrom, self.source, self.type, self.start, self.end,
                                                       str(self.score), self.strand, str(self.quality),
                                                       attribute_string)


class CollectionGFF(Collection):

    def read(self, input_file):
        self.metadata = MetadataGFF()
        self.records = []
        with open(input_file, "r") as in_fd:
            for line in in_fd:
                if line[0:2] == "##":
                    self.metadata.metadata.append(line[2:].strip())
                    continue
                elif line[0] == "#":
                    self.metadata.metadata.append(line[1:].strip())
                    continue
                tmp = line.strip().split("\t")
                #print(tmp)
                tmp[3] = int(tmp[3])
                tmp[4] = int(tmp[4])
                if tmp[5] != '.':
                    tmp[5] = float(tmp[5])
                if tmp[7] != '.':
                    tmp[7] = float(tmp[7])
                tmp[8] = tmp[8].split(";")
                tmp[8] = dict(map(lambda s: s.split("="), tmp[8]))

                self.records.append(RecordGFF(*tmp))

    @staticmethod
    def gff_simple_generator(input_file):

        with open(input_file, "r") as in_fd:
            for line in in_fd:
                if line[0] == "#":
                    continue
                tmp = line.strip().split("\t")
                #print(tmp)
                tmp[3] = int(tmp[3])
                tmp[4] = int(tmp[4])
                if tmp[5] != '.':
                    tmp[5] = float(tmp[5])
                if tmp[7] != '.':
                    tmp[7] = float(tmp[7])
                tmp[8] = tmp[8].split(";")
                tmp[8] = dict(map(lambda s: s.split("="), tmp[8]))
                yield tmp

    @staticmethod
    def gff_simple_record_generator(input_file):

        with open(input_file, "r") as in_fd:
            for line in in_fd:
                if line[0] == "#":
                    continue
                tmp = line.strip().split("\t")
                #print(tmp)
                tmp[3] = int(tmp[3])
                tmp[4] = int(tmp[4])
                if tmp[5] != '.':
                    tmp[5] = float(tmp[5])
                if tmp[7] != '.':
                    tmp[7] = float(tmp[7])
                tmp[8] = tmp[8].split(";")
                tmp[8] = dict(map(lambda s: s.split("="), tmp[8]))

                yield RecordGFF(*tmp)
