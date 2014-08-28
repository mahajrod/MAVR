#!/usr/bin/env python
from collections import Iterable
from Parser.Abstract import Record, Collection, Metadata
from Parser.VCF import CollectionVCF


class RecordCCF(Record, Iterable):

    def _check_chrom(self, vcf_records_list):
        #print (vcf_records_list)
        chrom = vcf_records_list[0].chrom
        for record in vcf_records_list:
            if record.chrom != chrom:
                raise ValueError("Records from different regions in same cluster")
        return chrom

    def __init__(self, id=None, chrom=None, size=None, start=None, end=None, description=None, flags=None,
                 vcf_records_list=None, bad_vcf_records=0, from_records=True, from_records_flag_mode="all"):
        # possible flags:
        # IP - indel(s) are present in record
        # BR - record is located in bad region
        self.records = vcf_records_list
        if from_records:
            self.chrom = self._check_chrom(vcf_records_list)
            self.size = len(vcf_records_list)
            self.start = vcf_records_list[0].pos
            self.end = vcf_records_list[-1].pos - 1 + \
                       max(map(lambda x: len(x), vcf_records_list[-1].alt_list + [vcf_records_list[-1].ref]))

            self.flags = set([])
            # possible from_records_flag_mode:
            # all - record to be counted as 'with flag' must have all flags from flags_list
            # one - record to be counted as 'with flag' must have at least one flag from flags_list
            if from_records_flag_mode == "one":
                for record in self:
                    self.flags |= record.flags
            elif from_records_flag_mode == "all":
                tmp_set = self.records[0].flags
                for record in self:
                    #print(tmp_set)
                    tmp_set &= record.flags
                self.flags |= tmp_set

            for record in self.records:
                if record.check_indel():
                    self.flags.add("IP")
                    break
            self.description = description

            if id:
                self.id = id
            else:
                self.id = "CL_%s_%i" % (self.chrom, self.start)
        else:
            if id:
                self.id = id
            else:
                self.id = "CL_%s_%i" % (self.chrom, self.start)
            self.chrom = chrom
            self.size = size
            self.start = start
            self.end = end
            self.description = description
            self.mean_dist = None
            self.flags = set(flags)
        self.len = self.end - self.start + 1
        self.bad_records = bad_vcf_records
        if bad_vcf_records > 0:
            self.flags.add("BR")

    def __len__(self):
        return self.len

    def __iter__(self):
        for record in self.records:
            yield record

    def __str__(self):
        attributes_string = "Size=%i;Bad_records=%i" % (self.size, self.bad_records)
        if self.flags:
            attributes_string += ";" + ";".join(self.flags)
        if self.description:
            attributes_string += ";" + ";".join(["%s=%s" % (key, str(self.description[key]))
                                                 for key in self.description])

        cluster_string = ">%s\t%s\t%i\t%i\t%s" % (self.id, self.chrom, self.start, self.end, attributes_string)
        return cluster_string + "\n\t" + "\n\t".join([str(record) for record in self.records])

    def check_location(self, bad_region_collection_gff):
        for bad_region in bad_region_collection_gff:
            for variant in self:
                if variant.chrom != bad_region.chrom:
                    continue
                if bad_region.start <= variant.pos <= bad_region.end:
                    self.bad_records += 1

        if self.bad_records > 0:
            self.flags.add("BR")

    def get_location(self, record_dict):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        if not self.description:
            self.description = {}
        if "Loc" not in self.description:
            self.description["Loc"] = set([])
        for variant in self:
            if "Loc" in variant.description:
                self.description["Loc"] |= set(variant.description["Loc"])
            for feature in record_dict[self.chrom].features:
                if (variant.pos - 1) in feature:
                    self.description["Loc"].add(feature.type)
                for sub_feature in feature.sub_features:
                    if (variant.pos - 1) in sub_feature:
                        self.description["Loc"].add(sub_feature.type)


class MetadataCCF(Metadata):

    def __init__(self, samples, metadata={}):
        self.samples = samples      #list
        self.metadata = metadata

    def __str__(self):
        metadata_string = "##Samples=" + ",".join(self.samples)
        if self.metadata:
            metadata_string += "\n##" + "\n##".join(["%s=%s" % (key, self.metadata[key]) for key in self.metadata])
        return metadata_string


class CollectionCCF(Collection):

    def read(self, input_file):
        # TODO: write read from ccf file; possible replace ccf by bed file
        pass

    def filter_by_size(self, min_size=3):
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            if record.size >= min_size:
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)
        return CollectionCCF(record_list=filtered_records), \
               CollectionCCF(record_list=filtered_out_records)

    def filter_by_flags(self, white_flag_list=[], black_flag_list=[]):
        filtered_records = []
        filtered_out_records = []
        white_list = set(white_flag_list)
        black_list = set(black_flag_list)
        for record in self.records:
            if white_list:
                if (white_list & record.flags) and not (black_list & record.flags):
                    filtered_records.append(record)
                else:
                    filtered_out_records.append(record)
            else:
                if black_list & record.flags:
                    filtered_out_records.append(record)
                else:
                    filtered_records.append(record)
        return CollectionCCF(record_list=filtered_records), \
               CollectionCCF(record_list=filtered_out_records)

    def check_record_location(self, bad_region_collection_gff):
        for record in self:
            record.check_location(bad_region_collection_gff)

    def get_collection_vcf(self, metadata, header):
        vcf_records = []
        for cluster in self:
            vcf_records += cluster.records
        return CollectionVCF(metadata=metadata, header_list=header, record_list=vcf_records)