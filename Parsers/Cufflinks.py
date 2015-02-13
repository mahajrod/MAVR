__author__ = 'mahajrod'
import re
from Parsers.Abstract import Record, Header, Metadata, Collection


class RecordFPKMTracing(Record):
    def __init__(self, chrom, pos, length, tracking_id, class_code, nearest_ref_id, gene_id, gene_short_name, tss_id,
                 coverage, FPKM, FPKM_conf_lo, FPKM_conf_hi, FPKM_status, description=None, flags=None):

        Record.__init__(self, chrom, pos,  description=description, flags=flags)
        # note that in .fpkm_tracking file locus coordinates are in python notation - starts from 0 and end nt is not included
        # probably, it is an error - check for every cufflinks version independently
        self.chrom = chrom                              # str
        self.length = length                            # int
        self.tracking_id = tracking_id                  # str
        self.class_code = class_code                    # str
        self.nearest_ref_id = nearest_ref_id            # str
        self.gene_id = gene_id                          # str
        self.gene_short_name = gene_short_name          # str
        self.tss_id = tss_id                            # str
        self.coverage = coverage                        # float or "-"
        self.FPKM = FPKM                                # float
        self.FPKM_conf_lo = FPKM_conf_lo                # float
        self.FPKM_conf_hi = FPKM_conf_hi                # float
        self.FPKM_status = FPKM_status                  # str

    def __str__(self):
        #tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	FPKM	FPKM_conf_lo	FPKM_conf_hi	FPKM_status
        return "\t".join([str(x) for x in (self.tracking_id, self.class_code, self.nearest_ref_id, self.gene_id,
                                           self.gene_short_name, self.tss_id, "%s:%i-%i" % (self.chrom, self.pos,
                                                                                            self.pos + self.length - 1),
                                           self.length, self.coverage, self.FPKM, self.FPKM_conf_lo, self.FPKM_conf_hi,
                                           self.FPKM_status)])


class HeaderFPKMTracking(Header):

    def __init__(self, header_string=None):
        if header_string is None:
            self.header = "tracking_id\tclass_code\tnearest_ref_id\tgene_id\tgene_short_name\ttss_id\tlocus\tlength\tcoverage\tFPKM\tFPKM_conf_lo\tFPKM_conf_hi\tFPKM_status"
        else:
            self.header = header_string

    def __str__(self):
        return self.header


class CollectionFPKMTracking(Collection):
    def __init__(self, metadata=None, record_list=None, header=None, input_file=None, from_file=False):
        Collection.__init__(self, metadata=metadata, record_list=record_list,
                            header=header if header else HeaderFPKMTracking(),
                            input_file=input_file, from_file=from_file)

    def read(self, input_file, input_in_python_coord_notation=True):
        self.metadata = None
        self.records = []
        with open(input_file, "r") as in_fd:
            self.header = in_fd.readline().strip()
            for line in in_fd:
                line_list = re.split("\t", line.strip())# line.strip().split("\t")
                if line[0] == "\t":
                    line_list.insert(0, "")

                chrom, pos = line_list[6].split(":")
                pos, end = [int(x) for x in pos.split("-")]

                if input_in_python_coord_notation:
                    pos += 1

                length = end - pos + 1
                cov_fpkm = [x if x == "-" else float(x) for x in line_list[8:-1]]
                packed_args = [chrom, pos, length] + line_list[:6] + cov_fpkm + [line_list[-1]]
                self.records.append(RecordFPKMTracing(*packed_args))

    def filter_by_expression(self, expression):
        filtered_records, filtered_out_records = self.filter_records_by_expression(expression)
        return CollectionFPKMTracking(metadata=self.metadata, header=self.header, record_list=filtered_records,
                                      from_file=False), \
               CollectionFPKMTracking(metadata=self.metadata, header=self.header, record_list=filtered_out_records,
                                      from_file=False)

if __name__ == "__main__":
    collection = CollectionFPKMTracking(from_file=True, input_file="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/Yeast_RNA_Seq/S1_Nagal_clout/genes.fpkm_tracking")
    collection.write("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/Yeast_RNA_Seq/S1_Nagal_clout/check.fpkm_tracking")