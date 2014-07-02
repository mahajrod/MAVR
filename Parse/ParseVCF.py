#!/usr/bin/python2

from collections import OrderedDict

import numpy as np


class RecordVCF():
    def __init__(self, chrom, pos, id, ref, alt_list, qual, filter_list, info_dict, samples_list):
        self.chrom = chrom                              #str
        self.pos = pos                                  #int
        self.id = id                                    #str
        self.ref = ref                                  #str
        self.alt_list = alt_list                        #list, entries are strings
        self.qual = qual                                #real or "."
        self.filter_list = sorted(filter_list)          #list, entries are strings
        self.info_dict = info_dict                      #dict
        #self.format_list = format_list                  #list, enties are strings
        self.samples_list = samples_list                #list entries are dicts with keys from format_list and
                                                        #values are lists

        #TODO: add data check
        #TODO: check parsing of files with several samples
    def __str__(self):
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        return self.string_form()

    def string_form(self):
        alt_string = ",".join(self.alt_list)
        filter_string = ";".join(self.filter_list)
        info_string = ";".join([key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key]))
                                for key in sorted(list(self.info_dict.keys()))])
        #format_string = ":".join(self.format_list)
        format_string = ":".join(self.samples_list[0].keys())
        samples_string = "\t".join([":".join([",".join(map(lambda x: str(x), sample[key])) for key in sample.keys()]) for sample in self.samples_list])

        #samples_string = "\t".join([":".join([",".join(str(sample[key])) for key in sample.keys()]) for sample in self.samples_list])

        return '\t'.join(map(lambda x: str(x), [self.chrom, self.pos, self.id, self.ref, alt_string,
                                                self.qual, filter_string, info_string, format_string, samples_string]))


class CollectionVCF():
    def __add_record(self, line):
        line_list = line.strip().split("\t")
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        position = int(line_list[1])
        quality = "."
        if quality != line_list[5]:
            quality = float(line_list[5])
        alt_list = line_list[4].split(",")
        filter_list = line_list[6].split(",")          #list, entries are strings

        info_tuple_list = [self.__split_by_equal_sign(entry) for entry in line_list[7].split(";")]
        info_dict = {}
        for entry in info_tuple_list:
            if self.metadata["INFO"][entry[0]]["Type"] == "Flag":
                info_dict[entry[0]] = []
            elif self.metadata["INFO"][entry[0]]["Type"] == "Integer":
                info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
            elif self.metadata["INFO"][entry[0]]["Type"] == "Float":
                info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
            else:
                info_dict[entry[0]] = entry[1].split(",")
        samples_list = []

        for sample_string in line_list[9:]:
            sample_dict = OrderedDict({})
            for key, value_list in zip(line_list[8].split(":"), sample_string.split(":")):
                #print (key, value_list)
                if self.metadata["FORMAT"][key]["Type"] == "Integer":
                    sample_dict[key] = list(map(lambda x: int(x), value_list.split(",")))
                elif self.metadata["FORMAT"][key]["Type"] == "Float":
                    sample_dict[key] = list(map(lambda x: float(x), value_list.split(",")))
                else:
                    sample_dict[key] = value_list.split(",")
            #print(sample_dict)
            samples_list.append(sample_dict)
        #print(samples_list)


        self.records.append(RecordVCF(line_list[0], position, line_list[2], line_list[3], alt_list, quality, filter_list,
                                      info_dict, samples_list))

        pass

    def __split_by_equal_sign(self, string):
        index = string.index("=")
        return string[:index], string[index+1:]

    def __split_by_comma_sign(self, string):
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == ",":
                index_list.append(i)
            i += 1
        index_list.append(len(string) - 1)
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def __add_metadata(self, line):
        key, value = self.__split_by_equal_sign(line[2:].strip())
        if value[0] == "<" and value[-1] == ">":
            #checking is value a list or no
            value = self.__split_by_comma_sign(value[1:-1])
            #print(value)
            #parse in suppose that first parameter in value list is ID
            value_id = self.__split_by_equal_sign(value[0])[1]
            #print(value[1:])
            value = dict(self.__split_by_equal_sign(entry) for entry in value[1:])
            #print(value_id, value)
            if key not in self.metadata:
                self.metadata[key] = {}
            self.metadata[key][value_id] = value
        else:
            self.metadata[key] = value

    def __init__(self, metadata=None, record_list=None, vcf_file=None, from_file=True):
        if from_file:
            self.metadata = OrderedDict({})
            self.records = []
            with open(vcf_file, "r") as fd:
                for line in fd:
                    if line[:2] != "##":
                        self.header_list = line[1:].strip().split("\t")
                        self.samples = self.header_list[8:]
                        break
                    #print(line)
                    self.__add_metadata(line)
                for line in fd:
                    self.__add_record(line)
        else:
            self.metadata = metadata
            self.records = record_list

    def write(self, output_file):
        pass

    def __len__(self):
        return len(self.records)

    def record_coordinates(self, black_list=[], white_list=[]):
        #return dictionary, where keys are seqids and values numpy arrays of SNV coordinates
        sequence_varcoord_dict = {}
        for record in self.records:
            if black_list and (record.chrom in black_list):
                continue
            if white_list and (record.chrom not in white_list):
                continue
            if record.chrom not in sequence_varcoord_dict:
                sequence_varcoord_dict[record.chrom] = [record.pos]
            else:
                sequence_varcoord_dict[record.chrom].append(record.pos)
        for chrom in sequence_varcoord_dict:
            sequence_varcoord_dict[chrom] = np.array(sequence_varcoord_dict[chrom])
        return sequence_varcoord_dict

    def stack_regions(self):
        #Not implemented yet
        #TODO: write
        pass
        """
        #assume that records inside regions are sorted by start coordinate
        sequence_varcoord_dict = self.record_coordinates()
        staked = np.array([])
        if "sequence-region" in self.metadata:
            shift_dict = {}
            shift = 0
            for seqid in sorted(list(self.metadata["sequence-region"].keys())):
                shift_dict[seqid] = shift
                sequence_varcoord_dict[seqid] += shift - self.metadata["sequence-region"][seqid][0] + 1
                staked = np.hstack((staked, sequence_varcoord_dict[seqid]))
                length = self.metadata["sequence-region"][seqid][1] - \
                         self.metadata["sequence-region"][seqid][0] + 1
                shift += length
        else:
            #TODO:implement stacking for files without sequence-region metadata
            pass
        return {"All": staked}, shift_dict
        """

    def rainfall_plot(self, output_file):
        #TODO: write rainfall plot
        pass

"""
class MutRecord(object):
    'docstring for MutRecord'
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
"""

if __name__ == "__main__":
    Mutations = CollectionVCF(vcf_file="/home/winstorage/old/GATK_vcf/filtered_snps/Sample_1_GATK_filtered_snps.vcf",
                              from_file=True)
    print(Mutations.metadata)
    print(len(Mutations))
    #print(Mutations.records[9])