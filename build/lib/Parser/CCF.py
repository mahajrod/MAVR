#!/usr/bin/env python

from collections import Iterable
import numpy as np
import matplotlib.pyplot as plt

from Parser.Abstract import Record, Collection, Metadata
from Parser.VCF import CollectionVCF


class RecordCCF(Record, Iterable):
    @staticmethod
    def _check_chrom(vcf_records_list):
        #print (vcf_records_list)
        #print(vcf_records_list.records)
        chrom = vcf_records_list[0].chrom
        for record in vcf_records_list:
            if record.chrom != chrom:
                raise ValueError("Records from different regions in same cluster")
        return chrom

    def __init__(self, id=None, chrom=None, size=None, start=None, end=None, description=None, flags=None,
                 collection_vcf=None, bad_vcf_records=0, from_records=True, subclusters=None, features=None,
                 ):
        # possible flags:
        # IP - indel(s) are present in record
        # BR - record is located in bad region
        self.records = collection_vcf
        self.subclusters = subclusters
        if from_records:
            self.chrom = self._check_chrom(collection_vcf)
            self.size = len(collection_vcf)
            self.start = collection_vcf[0].pos
            self.end = collection_vcf[-1].pos - 1 + \
                       max(map(lambda x: len(x), collection_vcf[-1].alt_list + [collection_vcf[-1].ref]))

            for record in self.records:
                if record.check_indel():
                    self.flags.add("IP")
                    break
            self.description = description if description else {}
            self.features = features
            if id:
                self.id = id
            else:
                self.id = "CL_%s_%i" % (self.chrom, self.start)
            self.flags = set(flags) if flags else set()
            self.bad_records = 0
            for record in self.records:
                if "BR" in record.flags:
                    self.bad_records += 1
        else:
            if id:
                self.id = id
            else:
                self.id = "CL_%s_%i" % (self.chrom, self.start)
            self.chrom = chrom
            self.size = size
            self.start = start
            self.end = end
            self.features = features
            self.description = description
            self.mean_dist = None
            self.flags = set(flags) if flags else set()
            self.bad_records = bad_vcf_records
        self.len = self.end - self.start + 1
        if bad_vcf_records > 0:
            self.flags.add("BR")
        else:
            self.flags - set(["BR"])

    def __len__(self):
        return self.size

    def __iter__(self):
        for record in self.records:
            yield record

    def __str__(self):
        attributes_string = "Size=%i;Bad_records=%i" % (self.size, self.bad_records)
        if self.flags:
            attributes_string += ";" + ";".join(self.flags)
        if self.description:
            attributes_string += ";" + ";".join(["%s=%s" % (key, ",".join(map(lambda x: str(x), self.description[key])))
                                                 for key in self.description])
        if self.subclusters != None:
            #print(self.subclusters)
            attributes_string += ";Subclusters=" + ",".join(map(lambda x: str(x), self.subclusters))

        cluster_string = ">%s\t%s\t%i\t%i\t%s" % (self.id, self.chrom, self.start, self.end, attributes_string)
        return cluster_string + "\nVariants\n\t" + "\n\t".join([str(record) for record in self.records])

    def reset_flags(self, flags_to_reset='all'):
        # flags_to_reset can take values 'all', set or list, or None
        # if all - all flags will be removed
        # if None - nothing will be done
        # if set/list - only flags in set/list will be removed
        if flags_to_reset == 'all':
            self.flags = set([])
        elif flags_to_reset:
            flags_to_remove = set(flags_to_reset)
            self.flags = self.flags - flags_to_remove

    def check_flags(self, flag_list, mismatch_list=[], expression_list=[], remove_mismatch_list=[],
                    flags_to_reset=None, mode="all", min_cluster_size=[]):
        # TODO: fix changes in  loc in description after removal
        # at moment subclustering resets after check if some record was removed

        # allow_flag_mismatch - number of allowed mismaches. if zero or None no mismatches are allowed
        self.reset_flags(flags_to_reset)
        if self.flags is None:
            self.flags = set([])


        mismatches = mismatch_list if mismatch_list else [0 for flag in flag_list]
        min_cluster = min_cluster_size if min_cluster_size else [3 for flag in flag_list]
        remove_list = remove_mismatch_list if remove_mismatch_list else[False for flag in flag_list]
        # possible from_records_flag_mode:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list

        # option flags_to_consider works only for remove_mismatch_records
        if mode == "one":
            for record in self:
                self.flags |= record.flags
        elif mode == "all":
            record_to_remove_dict = dict([(flag, set([])) for flag in flag_list])
            mismatch_count_dict = dict([(flag, 0) for flag in flag_list])
            index = 0
            for record in self:
                expressions = [eval(expression) for expression in expression_list] \
                    if expression_list else [True for flag in flag_list]
                for flag, expression, min_size in zip(flag_list, expressions, min_cluster):
                    if self.size < min_size:
                        continue
                    if flag not in record.flags:
                        if mismatch_count_dict[flag] is not None:
                            mismatch_count_dict[flag] += 1
                        if expression:
                            record_to_remove_dict[flag].add(index)
                        else:
                            # if expression is not True flag will not be set, None is used to indicate it
                            mismatch_count_dict[flag] = None
                index += 1
            for flag, mismatch in zip(flag_list, mismatches):
                #print(mismatch_count_dict, record_to_remove_dict)
                if (mismatch_count_dict[flag] is not None) and mismatch_count_dict[flag] <= mismatch:
                    self.flags.add(flag)
                    if mismatch_count_dict[flag] > 0:
                        self.description["N%sR" % flag] = [mismatch_count_dict[flag]]
                else:
                    record_to_remove_dict[flag].clear()
            records_to_remove = set([])
            for flag, remove in zip(flag_list, remove_list):
                if remove and (mismatch_count_dict[flag] is not None):
                    records_to_remove |= record_to_remove_dict[flag]
            records_to_remove = sorted(records_to_remove, reverse=True)
            #print(records_to_remove)
            if records_to_remove:
                #print(self)
                for index in records_to_remove:
                    self.records.pop(index)
                    self.subclusters = np.delete(self.subclusters, index)

                self.__init__(collection_vcf=self.records, description=self.description,
                              flags=self.flags, subclusters=self.subclusters,
                              from_records=True)

    def check_location(self, bad_region_collection_gff):
        self.bad_records = 0
        for variant in self:
            if "BR" in variant.flags:
                self.bad_records += 1
        if self.bad_records > 0:
            self.flags.add("BR")

    def get_location(self, record_dict, key="Loc"):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        if not self.description:
            self.description = {}
        if not self.features:
            self.features = []
        if key not in self.description:
            self.description[key] = set([])
        for variant in self:
            if key in variant.description:
                self.description[key] |= set(variant.description[key])
            for feature in record_dict[self.chrom].features:
                if (variant.pos - 1) in feature:
                    self.features.append(feature)
                    self.description[key].add(feature.type)
                for sub_feature in feature.sub_features:
                    if (variant.pos - 1) in sub_feature:
                        self.description[key].add(sub_feature.type)

    def subclustering(self,
                      method="inconsistent",
                      threshold=0.8,
                      cluster_distance='average'):
        tmp = self.records.get_clusters(extracting_method=method,
                                        threshold=threshold,
                                        cluster_distance=cluster_distance,
                                        split_by_regions=False,
                                        draw_dendrogramm=False,
                                        return_collection=False,
                                        write_inconsistent=False,
                                        write_correlation=False)
        self.subclusters = tmp[tmp.keys()[0]]

    def adjust(self, border_limit=None, min_size_to_adjust=2, remove_border_subclusters=False, remove_size_limit=1):
        # skip adjustment for clusters with 3 or less mutations
        if (self.size < min_size_to_adjust) or (self.subclusters is None):
            return -1
        limit = border_limit if border_limit else len(self.subclusters)
        for i in range(0, limit):
            if self.subclusters[i] == self.subclusters[0]:
                left_subcluster_end = i
            else:
                break
        # exit if cluster doesnt have subclusters
        if left_subcluster_end == len(self.subclusters) - 1:
            return 1

        for i in range(-1, -limit - 1, -1):
            if self.subclusters[i] == self.subclusters[-1]:
                right_subcluster_start = i
            else:
                break

        if remove_border_subclusters:
            start = left_subcluster_end + 1 if left_subcluster_end < remove_size_limit else 0
            end = right_subcluster_start if right_subcluster_start >= -remove_size_limit else len(self.subclusters)

            self.__init__(collection_vcf=CollectionVCF(record_list=self.records.records[start:end], from_file=False),
                          subclusters=self.subclusters[start:end], from_records=True)


class MetadataCCF(Metadata):

    def __init__(self, samples, vcf_metadata=None, vcf_header=None, metadata={}):
        self.samples = samples      #list
        self.metadata = metadata
        self.vcf_metadata = vcf_metadata
        self.vcf_header = vcf_header

    def __str__(self):
        metadata_string = None
        if self.vcf_metadata:
            metadata_string = "#VCF_METADATA START\n" + str(self.vcf_metadata) + "\n#VCF_METADATA END"
        if self.vcf_header:
            metadata_string = metadata_string + "\n#VCF_HEADER\n" + str(self.vcf_header) \
                if metadata_string else "#VCF_HEADER\n" + str(self.vcf_header)
        if self.metadata:
            metadata_string += "\n##" + "\n##".join(["%s=%s" % (key, self.metadata[key]) for key in self.metadata])
        return metadata_string


class CollectionCCF(Collection):

    def read(self, input_file):
        # TODO: write read from ccf file
        pass

    def filter_by_expression(self, expression):
        filtered_records, filtered_out_records = self.filter_records_by_expression(expression)
        return CollectionCCF(metadata=self.metadata, record_list=filtered_records,
                             from_file=False), \
               CollectionCCF(metadata=self.metadata, record_list=filtered_out_records,
                             from_file=False)

    def filter_by_size(self, min_size=3):
        return self.filter_by_expression("record.size >= %i" % min_size)

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
        return CollectionCCF(metadata=self.metadata, record_list=filtered_records), \
               CollectionCCF(metadata=self.metadata, record_list=filtered_out_records)

    def check_record_location(self, bad_region_collection_gff):
        for record in self:
            record.check_location(bad_region_collection_gff)

    def adjust(self, border_limit=None, min_size_to_adjust=2, remove_border_subclusters=False, remove_size_limit=1):
        for record in self:
            record.adjust(border_limit=border_limit, min_size_to_adjust=min_size_to_adjust,
                          remove_border_subclusters=remove_border_subclusters, remove_size_limit=remove_size_limit)

    def subclustering(self,
                      method="inconsistent",
                      threshold=0.8,
                      cluster_distance='average'):
        for record in self:
            if len(record) < 3:
                continue
            #print(record)
            record.subclustering(method=method,
                                 threshold=threshold,
                                 cluster_distance=cluster_distance)
    """
    def get_collection_vcf(self, metadata, header):
        vcf_records = []
        for cluster in self:
            vcf_records += cluster.records
        return CollectionVCF(metadata=metadata, header=header, record_list=vcf_records)
    """
    def count(self):
        sizes = []
        for record in self:
            sizes.append(record.size)
        return sizes

    def statistics(self, filename="cluster_size_distribution.svg", title="Distribution of sizes of clusters",
                   dpi=150, figsize=(10, 10), facecolor="green"):
        plt.figure(1, dpi=dpi, figsize=figsize)
        plt.subplot(1, 1, 1)
        plt.suptitle(title)
        counts = self.count()
        maximum = max(counts)
        bins = np.linspace(0, maximum, maximum)
        plt.hist(counts, bins, facecolor=facecolor)
        plt.xticks(np.arange(0, maximum, 1.0))
        plt.savefig(filename)
        plt.close()

    def check_flags(self, flag_list, mismatch_list=[], expression_list=[], remove_mismatch_list=[],
                    flags_to_reset=None, mode="all", min_cluster_size=[]):
        for record in self:
            record.check_flags(flag_list, mismatch_list=mismatch_list, expression_list=expression_list,
                               remove_mismatch_list=remove_mismatch_list, flags_to_reset=flags_to_reset, mode=mode,
                               min_cluster_size=min_cluster_size)

    def extract_vcf(self):
        vcf = CollectionVCF(metadata=self.metadata.vcf_metadata, record_list=[], header=self.metadata.vcf_header,
                            samples=self.metadata.samples, from_file=False)
        for record in self:
            """
            print(record)
            print(type(record))
            print(record.records)
            print(type(record.records))
            """
            vcf = vcf + record.records
        return vcf