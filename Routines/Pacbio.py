#!/usr/bin/env python

from Routines.Alignment import AlignmentRoutines


class PacBioRoutines(AlignmentRoutines):
    def __init__(self):
        AlignmentRoutines.__init__(self)

    def extract_not_collapsed_subreads_from_bam(self, subread_bam, ccs_bam, output_prefix):
        """
        ccs_read_names_file = "%s.ccs.reads.ids" % output_prefix
        ccs_read_names_prefix_file = "%s.ccs.reads.prefix.ids" % output_prefix

        SamtoolsV1.get_read_names(ccs_bam, ccs_read_names_file)

        sed_string = "sed \'s/ccs$//\' %s > %s " % (ccs_read_names_file, ccs_read_names_prefix_file)
        os.system(sed_string)
        read_name_list = IdList(filename=ccs_read_names_prefix_file)

        SamtoolsV1.get_reads_by_name(read_name_list, input_sam_fd, output_sam_fd, mode="include", search_mode="exact"):

        """


