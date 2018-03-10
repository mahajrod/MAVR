#!/usr/bin/env python

import re

from CustomCollections.GeneralCollections import IdList

from Parsers.Abstract import Collection

from Routines import AnnotationsRoutines


class PrimerEntryPrimer3:
    def __init__(self, penalty, seq, start, length, melting_temperature, gc_content,
                 self_any_th, self_end_th, hairpin_th, end_stability):
        self.penalty = penalty                              #float
        self.seq = seq                                      #str
        self.start = start                                  #int
        self.length = length                                #int
        self.melting_temperature = melting_temperature      #float
        self.gc_content = gc_content                        #float
        self.self_any_th = self_any_th                      #float
        self.self_end_th = self_end_th                      #float
        self.hairpin_th = hairpin_th                        #float
        self.end_stability = end_stability                  #float


class PrimerPairEntryPrimer3:
    def __init__(self, pair_id=None, product_size=None, penalty=None, compl_any_th=None, compl_end_th=None,
                 left_primer=None, right_primer=None):
        self.id = pair_id                       #int
        self.product_size = product_size        #int
        self.penalty = penalty                  #float
        self.compl_any_th = compl_any_th        #float
        self.compl_end_th = compl_end_th        #float
        self.left_primer = left_primer          #PrimerEntryPrimer3
        self.right_primer = right_primer        #PrimerEntryPrimer3

    def primer_string(self):
        string = ""

        string += "PRIMER_LEFT_%i_PENALTY=%f\n" % (self.id, self.left_primer.penalty)
        string += "PRIMER_RIGHT_%i_PENALTY=%f\n" % (self.id, self.right_primer.penalty)
        string += "PRIMER_LEFT_%i_SEQUENCE=%s\n" % (self.id, self.left_primer.seq)
        string += "PRIMER_RIGHT_%i_SEQUENCE=%s\n" % (self.id, self.right_primer.seq)
        string += "PRIMER_LEFT_%i=%i,%i\n" % (self.id, self.left_primer.start, self.left_primer.length)
        string += "PRIMER_RIGHT_%i=%i,%i\n" % (self.id, self.right_primer.start, self.right_primer.length)
        string += "PRIMER_LEFT_%i_TM=%f\n" % (self.id, self.left_primer.melting_temperature)
        string += "PRIMER_RIGHT_%i_TM=%f\n" % (self.id, self.right_primer.melting_temperature)
        string += "PRIMER_LEFT_%i_GC_PERCENT=%f\n" % (self.id, self.left_primer.gc_content)
        string += "PRIMER_RIGHT_%i_GC_PERCENT=%f\n" % (self.id, self.right_primer.gc_content)
        string += "PRIMER_LEFT_%i_SELF_ANY_TH=%f\n" % (self.id, self.left_primer.self_any_th)
        string += "PRIMER_RIGHT_%i_SELF_ANY_TH=%f\n" % (self.id, self.right_primer.self_any_th)
        string += "PRIMER_LEFT_%i_SELF_END_TH=%f\n" % (self.id, self.left_primer.self_end_th)
        string += "PRIMER_RIGHT_%i_SELF_END_TH=%f\n" % (self.id, self.right_primer.self_end_th)
        string += "PRIMER_LEFT_%i_HAIRPIN_TH=%f\n" % (self.id, self.left_primer.hairpin_th)
        string += "PRIMER_RIGHT_%i_HAIRPIN_TH=%f\n" % (self.id, self.right_primer.hairpin_th)
        string += "PRIMER_LEFT_%i_END_STABILITY=%f\n" % (self.id, self.left_primer.end_stability)
        string += "PRIMER_RIGHT_%i_END_STABILITY=%f\n" % (self.id, self.right_primer.end_stability)

        return string

    def __str__(self):
        if self.id is None:
            return None
        string = "PRIMER_PAIR_%i_PENALTY=%f\n" % (self.id, self.penalty)

        string += self.primer_string()

        string += "PRIMER_PAIR_%i_COMPL_ANY_TH=%f\n" % (self.id, self.compl_any_th) if self.compl_any_th else ""
        string += "PRIMER_PAIR_%i_COMPL_END_TH=%f\n" % (self.id, self.compl_end_th) if self.compl_end_th else ""
        string += "PRIMER_PAIR_%i_PRODUCT_SIZE=%i\n" % (self.id, self.product_size)

        return string

    def pcr_product_seq(self, seq):
        if self.left_primer and self.right_primer:
            return seq[self.left_primer.start:self.right_primer.start+1]
        return None

    def gaps_inside_pcr_product(self, seq, min_gap_len=5):
        if self.left_primer and self.right_primer:
            pcr_seq = self.pcr_product_seq(seq)
            gap_reg_exp = re.compile("N{%i,}" % min_gap_len, re.IGNORECASE)
            gaps_list = list(gap_reg_exp.finditer(pcr_seq))
            if gaps_list:
                return True
            return False
        return None


class RecordPrimer3:
    def __init__(self, seq_id, seq, target_start, target_len,
                 pick_left_primer, pick_internal_oligo, pick_right_primer, product_size_range,
                 left_primer_count, internal_oligo_count, right_primer_count, primer_pair_count, primer_pair_list=[],
                 left_primer_choice_description=None, right_primer_choice_description=None,
                 internal_oligo_choice_description=None, pair_choice_description=None, chrom=None,
                 chrom_pos_start=None, chrom_pos_end=None, monomer=None, total_monomer_number=None):

        self.id = seq_id                                                                    # str
        self.seq = seq                                                                      # str
        self.target_start = target_start                                                    # int
        self.target_len = target_len                                                        # int
        self.pick_left_primer = pick_left_primer                                            # bool
        self.pick_internal_oligo = pick_internal_oligo                                      # bool
        self.pick_right_primer = pick_right_primer                                          # bool
        self.product_size_range = product_size_range                                        # (int, int)
        self.left_primer_count = left_primer_count                                          # int
        self.internal_oligo_count = internal_oligo_count                                    # int
        self.right_primer_count = right_primer_count                                        # int
        self.primer_pair_count = primer_pair_count                                          # int
        self.left_primer_choice_description = left_primer_choice_description                # str
        self.right_primer_choice_description = right_primer_choice_description              # str
        self.internal_oligo_choice_description = internal_oligo_choice_description          # str
        self.pair_choice_description = pair_choice_description                              # str
        self.primer_pair_list = primer_pair_list                                            # list of PrimerPairEntryPrimer3

        self.chrom = chrom                                                                  # str
        self.chrom_pos_start = chrom_pos_start                                              # int
        self.chrom_pos_end = chrom_pos_end                                                  # int

        self.monomer = monomer                                                              # str
        self.total_monomer_number = total_monomer_number                                    # float

    def __str__(self):
        string = ""
        string += "SEQUENCE_ID=%s\n" % self.id
        string += "SEQUENCE_TEMPLATE=%s\n" % self.seq
        string += "SEQUENCE_TARGET=%i,%i\n" % (self.target_start, self.target_len)
        string += "PRIMER_PICK_LEFT_PRIMER=%i\n" % (1 if self.pick_left_primer else 0) if not (self.pick_left_primer is None) else ""
        string += "PRIMER_PICK_INTERNAL_OLIGO=%i\n" % (1 if self.pick_internal_oligo else 0) if not (self.pick_internal_oligo is None) else ""
        string += "PRIMER_PICK_RIGHT_PRIMER=%i\n" % (1 if self.pick_right_primer else 0) if not (self.pick_right_primer is None) else ""
        string += "PRIMER_PRODUCT_SIZE_RANGE=%i-%i\n" % (self.product_size_range[0], self.product_size_range[1])
        string += "PRIMER_INTERNAL_OLIGO_EXPLAIN=%s\n" % self.left_primer_choice_description if not (self.left_primer_choice_description is None) else ""
        string += "PRIMER_PAIR_EXPLAIN=%s\n" % self.internal_oligo_choice_description if not (self.internal_oligo_choice_description is None) else ""
        string += "PRIMER_RIGHT_EXPLAIN=%s\n" % self.right_primer_choice_description if not (self.right_primer_choice_description is None) else ""
        string += "PRIMER_PAIR_EXPLAIN=%s\n" % self.pair_choice_description if not (self.pair_choice_description is None) else ""
        string += "PRIMER_LEFT_NUM_RETURNED=%i\n" % self.left_primer_count
        string += "PRIMER_RIGHT_NUM_RETURNED=%i\n" % self.right_primer_count
        string += "PRIMER_INTERNAL_NUM_RETURNED=%i\n" % self.internal_oligo_count
        string += "PRIMER_PAIR_NUM_RETURNED=%i\n" % self.primer_pair_count

        for primer_pair in self.primer_pair_list:
            string += str(primer_pair)

        string += "=\n"

        return string

    def remove_primers_with_gaps_in_pcr_product(self, min_gap_len=5):
        bad_primers_index_list = []

        good_pairs_list = []
        for primer_pair_index in range(0, len(self.primer_pair_list)):
            if self.primer_pair_list[primer_pair_index].gaps_inside_pcr_product(self.seq, min_gap_len=min_gap_len):
                bad_primers_index_list.append(primer_pair_index)
            else:
                good_pairs_list.append(self.primer_pair_list[primer_pair_index])

        ok_pairs = len(good_pairs_list)

        if bad_primers_index_list:
            if self.pair_choice_description:
                self.pair_choice_description += ", gaps in pcr product %i, ok after filtering %i" % (len(bad_primers_index_list), ok_pairs)
            self.primer_pair_list = []
            for good_pair_index in range(0, len(good_pairs_list)):
                good_pairs_list[good_pair_index].id = good_pair_index
            self.primer_pair_list = good_pairs_list
            self.primer_pair_count = ok_pairs

        return len(bad_primers_index_list)

    def table_form(self):
        string = ""
        metadata_location_str = "#Location"
        if self.primer_pair_list:
            metadata_str = "#SeqeunceID\t%s" % self.id
            if self.chrom:
                metadata_location_str += "\t%s" % self.chrom
                if self.chrom_pos_start and self.chrom_pos_end:
                    metadata_location_str += ":%i-%i" % (self.chrom_pos_start, self.chrom_pos_end)
            metadata_seq_str = "#Sequence\t%s" % self.seq
            header_str = "Primer_pair"
            pcr_product_size_str = "PCR product size"

            primer_pair_penalty_str = "Primer pair penalty"
            left_primer_penalty_str = "Left primer penalty"
            right_primer_penalty_str = "Right primer penalty"

            left_primer_seq_str = "Left primer seq"
            right_primer_seq_str = "Right primer seq"
            left_primer_tm_str = "Left primer melting temp"
            right_primer_tm_str = "Right primer melting temp"
            left_primer_gc_str = "Left primer GC,%"
            right_primer_gc_str = "Right primer GC,%"

            left_primer_start_str = "Left primer start(1-based)"
            right_primer_start_str = "Right primer start(1-based)"

            left_primer_len_str = "Left primer length"
            right_primer_len_str = "Right primer length"

            for primer_pair in self.primer_pair_list:
                header_str += "\t%s" % str(primer_pair.id)
                pcr_product_size_str += "\t%i" % primer_pair.product_size

                primer_pair_penalty_str += "\t%f" % primer_pair.penalty
                left_primer_penalty_str += "\t%s" % primer_pair.left_primer.penalty
                right_primer_penalty_str += "\t%s" % primer_pair.right_primer.penalty

                left_primer_seq_str += "\t%s" % primer_pair.left_primer.seq
                right_primer_seq_str += "\t%s" % primer_pair.right_primer.seq

                left_primer_tm_str += "\t%f" % primer_pair.left_primer.melting_temperature
                right_primer_tm_str += "\t%f" % primer_pair.right_primer.melting_temperature

                left_primer_gc_str += "\t%f" % primer_pair.left_primer.gc_content
                right_primer_gc_str += "\t%f" % primer_pair.right_primer.gc_content

                left_primer_start_str += "\t%i" % (primer_pair.left_primer.start + 1)
                right_primer_start_str += "\t%i" % (primer_pair.right_primer.start + 1)

                left_primer_len_str += "\t%i" % primer_pair.left_primer.length
                right_primer_len_str += "\t%i" % primer_pair.right_primer.length

            for table_string in (metadata_str, metadata_location_str, metadata_seq_str,
                                 header_str, pcr_product_size_str, primer_pair_penalty_str,
                                 left_primer_penalty_str, right_primer_penalty_str,
                                 left_primer_seq_str, right_primer_seq_str,
                                 left_primer_tm_str, right_primer_tm_str,
                                 left_primer_gc_str, right_primer_gc_str,
                                 left_primer_len_str, right_primer_len_str,
                                 left_primer_start_str, right_primer_start_str):

                string += table_string + "\n"

        return string

    def table_form2_short(self):
        """
        header = "#STR_ID\tScaffold\tStart\tStop\tMonomer\tCopy_number\tPRIMER_PAIR_NUMBER\tPRIMER_PAIR_ID\tPCR_PRODUCT_SIZE\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\t" \
                 "PRIMER_LEFT_GC,%\tPRIMER_RIGHT_GC,%\tPRIMER_LEFT_SEQ\tPRIMER_RIGHT_SEQ\tTARGET_SEQ\n"
        """
        table_form2_short = ""
        for primer_pair in self.primer_pair_list:

            table_form2_short += "%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\n" % (self.id,
                                                                                                    self.chrom,
                                                                                                    self.chrom_pos_start,
                                                                                                    self.chrom_pos_end,
                                                                                                    str(self.monomer) if not(self.monomer is None) else ".",
                                                                                                    str(self.total_monomer_number) if not(self.total_monomer_number is None) else ".",
                                                                                                    self.primer_pair_count,
                                                                                                    primer_pair.id,
                                                                                                    primer_pair.product_size,
                                                                                                    primer_pair.left_primer.melting_temperature,
                                                                                                    primer_pair.right_primer.melting_temperature,
                                                                                                    primer_pair.left_primer.gc_content,
                                                                                                    primer_pair.right_primer.gc_content,
                                                                                                    primer_pair.left_primer.seq,
                                                                                                    primer_pair.right_primer.seq)
        return table_form2_short

    def table_form2(self):
        """
        header = "#STR_ID\tScaffold\tStart\tStop\tMonomer\tCopy_number\tPRIMER_PAIR_NUMBER\tPRIMER_PAIR_ID\tPCR_PRODUCT_SIZE\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\t" \
                 "PRIMER_LEFT_GC,%\tPRIMER_RIGHT_GC,%\tPRIMER_LEFT_SEQ\tPRIMER_RIGHT_SEQ\tTARGET_SEQ\t" \
                 "PRIMER_LEFT_SELF_ANY_TH\tPRIMER_RIGHT_SELF_ANY_TH\tPRIMER_LEFT_SELF_END_TH\tPRIMER_RIGHT_SELF_END_TH\t" \
                 "PRIMER_LEFT_HAIRPIN_TH\tPRIMER_RIGHT_HAIRPIN_TH\tPRIMER_LEFT_END_STABILITY\tPRIMER_RIGHT_END_STABILITY\tREGION_SEQ\n"


        """
        table_form2 = ""
        for primer_pair in self.primer_pair_list:

            table_form2 += "%s\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\t%.3f\t%.3f\t" \
                           "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n" % (self.id,
                                                                         self.chrom,
                                                                         self.chrom_pos_start,
                                                                         self.chrom_pos_end,
                                                                         str(self.monomer) if not(self.monomer is None) else ".",
                                                                         str(self.total_monomer_number) if not(self.total_monomer_number is None) else ".",
                                                                         self.primer_pair_count,
                                                                         primer_pair.id,
                                                                         primer_pair.product_size,
                                                                         primer_pair.left_primer.melting_temperature,
                                                                         primer_pair.right_primer.melting_temperature,
                                                                         primer_pair.left_primer.gc_content,
                                                                         primer_pair.right_primer.gc_content,
                                                                         primer_pair.left_primer.seq,
                                                                         primer_pair.right_primer.seq,
                                                                         primer_pair.left_primer.self_any_th,
                                                                         primer_pair.right_primer.self_any_th,
                                                                         primer_pair.left_primer.self_end_th,
                                                                         primer_pair.right_primer.self_end_th,
                                                                         primer_pair.left_primer.hairpin_th,
                                                                         primer_pair.right_primer.hairpin_th,
                                                                         primer_pair.left_primer.end_stability,
                                                                         primer_pair.right_primer.end_stability,
                                                                         self.seq)

        return table_form2

    def alignments_string(self, segment_length=120, left_primer_symbol=">",
                          target_symbol="*", right_primer_symbol="<"):
        string = ""
        string += "#SeqeunceID\t%s\n" % self.id
        string += "#Location"
        if self.chrom:
            string += "\t%s" % self.chrom
            if self.chrom_pos_start and self.chrom_pos_end:
                string += ":%i-%i" % (self.chrom_pos_start, self.chrom_pos_end)
        string += "\n"
        #string += "#Sequence\t%s\n" % self.seq

        for primer_pair in self.primer_pair_list:
            string += "#Primer pair %i\n" % primer_pair.id
            string += "\n"

            location_list = [(primer_pair.left_primer.start,
                              primer_pair.left_primer.start + primer_pair.left_primer.length),

                             (self.target_start, self.target_start + self.target_len),

                             (primer_pair.right_primer.start - primer_pair.right_primer.length + 1,
                              primer_pair.right_primer.start + 1)]
            #print self.id
            #print location_list
            string += AnnotationsRoutines.draw_string_regions(self.seq, location_list,
                                                           [left_primer_symbol, target_symbol, right_primer_symbol],
                                                           overlap_symbol="#", line_per_record=False,
                                                           segment_length=segment_length,
                                                           num_of_spaces=3, num_of_space_lines=1, empty_symbol=" ")

        return string

    def table_form_and_alignments(self, segment_length=120, left_primer_symbol=">",
                          target_symbol="*", right_primer_symbol="<"):

        string = self.table_form()
        string += self.alignments_string(segment_length=segment_length, left_primer_symbol=left_primer_symbol,
                                         target_symbol=target_symbol, right_primer_symbol=right_primer_symbol)

        return string


class CollectionPrimer3(Collection):

    def __init__(self, record_list=None, primer3_file=None, from_file=True, id_based_location_dict=None, 
                 repeat_gff_file=None, id_description_entry="ID"):
        self.general_entry_list = ["SEQUENCE_ID",
                                   "SEQUENCE_TEMPLATE",
                                   "SEQUENCE_TARGET",
                                   "PRIMER_PICK_LEFT_PRIMER",
                                   "PRIMER_PICK_INTERNAL_OLIGO",
                                   "PRIMER_PICK_RIGHT_PRIMER",
                                   "PRIMER_PRODUCT_SIZE_RANGE",
                                   "PRIMER_LEFT_EXPLAIN",
                                   "PRIMER_RIGHT_EXPLAIN",
                                   "PRIMER_PAIR_EXPLAIN",
                                   "PRIMER_LEFT_NUM_RETURNED",
                                   "PRIMER_RIGHT_NUM_RETURNED",
                                   "PRIMER_INTERNAL_NUM_RETURNED",
                                   "PRIMER_PAIR_NUM_RETURNED"]

        self.primer_entry_prefix_list = ["PRIMER_LEFT", "PRIMER_RIGHT"]
        self.primer_entry_suffix_list = ["PENALTY",
                                         "SEQUENCE",
                                         "TM",
                                         "GC_PERCENT",
                                         "SELF_ANY_TH",
                                         "SELF_END_TH",
                                         "HAIRPIN_TH",
                                         "END_STABILITY"]

        self.primer_pair_prefix_list = ["PRIMER_PAIR"]
        self.primer_pair_suffix_list = ["PENALTY",
                                        "COMPL_ANY_TH",
                                        "COMPL_END_TH",
                                        "PRODUCT_SIZE"]

        id_based_dict = id_based_location_dict

        if repeat_gff_file:
            id_based_dict = AnnotationsRoutines.get_id_based_dict_from_gff(repeat_gff_file,
                                                                           id_entry=id_description_entry)

        if from_file:
            self.records = []
            with open(primer3_file, "r") as in_fd:
                for line in in_fd:
                    entry_dict = {}
                    lineeee = line
                    while lineeee[0] != "=":
                        line_list = lineeee.strip().split("=")
                        entry_dict[line_list[0]] = line_list[1]
                        lineeee = in_fd.next()

                    self._add_record(entry_dict, id_based_location_dict=id_based_dict)

        else:
            self.records = record_list

    def _add_record(self, entry_dict, id_based_location_dict=None):

        target_start, target_len = map(int, entry_dict["SEQUENCE_TARGET"].split(","))

        primer_pair_number = int(entry_dict["PRIMER_PAIR_NUM_RETURNED"]) if "PRIMER_PAIR_NUM_RETURNED" in entry_dict else None

        primer_pair_list = []
        
        if id_based_location_dict:
            #print id_based_location_dict
            chrom = id_based_location_dict[entry_dict["SEQUENCE_ID"]][0]
            start = int(id_based_location_dict[entry_dict["SEQUENCE_ID"]][1])
            end = int(id_based_location_dict[entry_dict["SEQUENCE_ID"]][2])

            monomer = id_based_location_dict[entry_dict["SEQUENCE_ID"]][3]["Pattern"][0] if "Pattern" in id_based_location_dict[entry_dict["SEQUENCE_ID"]][3] else None
            total_monomer_number = float(id_based_location_dict[entry_dict["SEQUENCE_ID"]][3]["N_copies"][0]) if "N_copies" in id_based_location_dict[entry_dict["SEQUENCE_ID"]][3] else None
        else:
            chrom = None
            start = None
            end = None
            monomer = None
            total_monomer_number = None
        
        for primer_pair_index in range(0, primer_pair_number):
            left_primer_start, left_primer_len = map(int, entry_dict["PRIMER_LEFT_%i" % primer_pair_index].split(","))
            right_primer_start, right_primer_len = map(int, entry_dict["PRIMER_RIGHT_%i" % primer_pair_index].split(","))

            left_primer = PrimerEntryPrimer3(float(entry_dict["PRIMER_LEFT_%i_PENALTY" % primer_pair_index]),
                                             entry_dict["PRIMER_LEFT_%i_SEQUENCE" % primer_pair_index],
                                             left_primer_start,
                                             left_primer_len,
                                             float(entry_dict["PRIMER_LEFT_%i_TM" % primer_pair_index]),
                                             float(entry_dict["PRIMER_LEFT_%i_GC_PERCENT" % primer_pair_index]),
                                             float(entry_dict["PRIMER_LEFT_%i_SELF_ANY_TH" % primer_pair_index]),
                                             float(entry_dict["PRIMER_LEFT_%i_SELF_END_TH" % primer_pair_index]),
                                             float(entry_dict["PRIMER_LEFT_%i_HAIRPIN_TH" % primer_pair_index]),
                                             float(entry_dict["PRIMER_LEFT_%i_END_STABILITY" % primer_pair_index]))

            right_primer = PrimerEntryPrimer3(float(entry_dict["PRIMER_RIGHT_%i_PENALTY" % primer_pair_index]),
                                              entry_dict["PRIMER_RIGHT_%i_SEQUENCE" % primer_pair_index],
                                              right_primer_start,
                                              right_primer_len,
                                              float(entry_dict["PRIMER_RIGHT_%i_TM" % primer_pair_index]),
                                              float(entry_dict["PRIMER_RIGHT_%i_GC_PERCENT" % primer_pair_index]),
                                              float(entry_dict["PRIMER_RIGHT_%i_SELF_ANY_TH" % primer_pair_index]),
                                              float(entry_dict["PRIMER_RIGHT_%i_SELF_END_TH" % primer_pair_index]),
                                              float(entry_dict["PRIMER_RIGHT_%i_HAIRPIN_TH" % primer_pair_index]),
                                              float(entry_dict["PRIMER_RIGHT_%i_END_STABILITY" % primer_pair_index]))

            primer_pair_list.append(PrimerPairEntryPrimer3(pair_id=primer_pair_index,
                                                           product_size=int(entry_dict["PRIMER_PAIR_%i_PRODUCT_SIZE" % primer_pair_index]),
                                                           penalty=float(entry_dict["PRIMER_PAIR_%i_PENALTY" % primer_pair_index]),
                                                           compl_any_th=float(entry_dict["PRIMER_PAIR_%i_COMPL_ANY_TH" % primer_pair_index]) if "PRIMER_PAIR_%i_COMPL_ANY_TH" in entry_dict else None,
                                                           compl_end_th=float(entry_dict["PRIMER_PAIR_%i_COMPL_END_TH" % primer_pair_index]) if "PRIMER_PAIR_%i_COMPL_END_TH" in entry_dict else None,
                                                           left_primer=left_primer,
                                                           right_primer=right_primer))

        record = RecordPrimer3(entry_dict["SEQUENCE_ID"],
                               entry_dict["SEQUENCE_TEMPLATE"],
                               target_start, target_len,
                               int(entry_dict["PRIMER_PICK_LEFT_PRIMER"]) if "PRIMER_PICK_LEFT_PRIMER" in entry_dict else None,
                               int(entry_dict["PRIMER_PICK_INTERNAL_OLIGO"]) if "PRIMER_PICK_INTERNAL_OLIGO" in entry_dict else None,
                               int(entry_dict["PRIMER_PICK_RIGHT_PRIMER"]) if "PRIMER_PICK_RIGHT_PRIMER" in entry_dict else None,
                               map(int, entry_dict["PRIMER_PRODUCT_SIZE_RANGE"].split("-")),
                               int(entry_dict["PRIMER_LEFT_NUM_RETURNED"]) if "PRIMER_LEFT_NUM_RETURNED" in entry_dict else None,
                               int(entry_dict["PRIMER_INTERNAL_NUM_RETURNED"]) if "PRIMER_INTERNAL_NUM_RETURNED" in entry_dict else None,
                               int(entry_dict["PRIMER_RIGHT_NUM_RETURNED"]) if "PRIMER_RIGHT_NUM_RETURNED" in entry_dict else None,
                               primer_pair_number,
                               primer_pair_list=primer_pair_list,
                               left_primer_choice_description=entry_dict["PRIMER_LEFT_EXPLAIN"] if "PRIMER_LEFT_EXPLAIN" in entry_dict else None,
                               right_primer_choice_description=entry_dict["PRIMER_RIGHT_EXPLAIN"] if "PRIMER_RIGHT_EXPLAIN" in entry_dict else None,
                               internal_oligo_choice_description=entry_dict["PRIMER_INTERNAL_OLIGO_EXPLAIN"] if "PRIMER_INTERNAL_OLIGO_EXPLAIN" in entry_dict else None,
                               pair_choice_description=entry_dict["PRIMER_PAIR_EXPLAIN"] if "PRIMER_PAIR_EXPLAIN" in entry_dict else None,
                               chrom=chrom,
                               chrom_pos_start=start,
                               chrom_pos_end=end,
                               monomer=monomer,
                               total_monomer_number=total_monomer_number)

        self.records.append(record)

    def write(self, out_file):
        with open(out_file, "w") as out_fd:
            for record in self.records:
                out_fd.write(str(record))

    def write_table_form(self, out_file):
        with open(out_file, "w") as out_fd:
            for record in self.records:
                out_fd.write(record.table_form())

    def write_alignments(self, out_file, segment_length=120, left_primer_symbol=">",
                         target_symbol="*", right_primer_symbol="<"):
        with open(out_file, "w") as out_fd:
            for record in self.records:
                out_fd.write(record.alignments_string(segment_length=segment_length,
                                                      left_primer_symbol=left_primer_symbol,
                                                      target_symbol=target_symbol,
                                                      right_primer_symbol=right_primer_symbol))

    def write_table_form_with_alignments(self, out_file, segment_length=120, left_primer_symbol=">",
                                         target_symbol="*", right_primer_symbol="<"):
        with open(out_file, "w") as out_fd:
            for record in self.records:
                out_fd.write(record.table_form_and_alignments(segment_length=segment_length,
                                                              left_primer_symbol=left_primer_symbol,
                                                              target_symbol=target_symbol,
                                                              right_primer_symbol=right_primer_symbol))

    def write_table_form2_short(self, out_file):

        header = "#STR_ID\tScaffold\tStart\tStop\tMonomer\tCopy_number\tPRIMER_PAIR_NUMBER\tPRIMER_PAIR_ID" \
                 "\tPCR_PRODUCT_SIZE\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\t" \
                 "PRIMER_LEFT_GC,%\tPRIMER_RIGHT_GC,%\tPRIMER_LEFT_SEQ\tPRIMER_RIGHT_SEQ\tTARGET_SEQ\n"

        with open(out_file, "w") as out_fd:
            out_fd.write(header)
            for record in self.records:
                out_fd.write(record.table_form2_short())

    def write_table_form2(self, out_file):

        header = "#STR_ID\tScaffold\tStart\tStop\tMonomer\tCopy_number\tPRIMER_PAIR_NUMBER\tPRIMER_PAIR_ID" \
                 "\tPCR_PRODUCT_SIZE\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\t" \
                 "PRIMER_LEFT_GC,%\tPRIMER_RIGHT_GC,%\tPRIMER_LEFT_SEQ\tPRIMER_RIGHT_SEQ\t" \
                 "PRIMER_LEFT_SELF_ANY_TH\tPRIMER_RIGHT_SELF_ANY_TH\tPRIMER_LEFT_SELF_END_TH\tPRIMER_RIGHT_SELF_END_TH\t" \
                 "PRIMER_LEFT_HAIRPIN_TH\tPRIMER_RIGHT_HAIRPIN_TH\tPRIMER_LEFT_END_STABILITY\tPRIMER_RIGHT_END_STABILITY\tTARGET_SEQ\n"

        with open(out_file, "w") as out_fd:
            out_fd.write(header)
            for record in self.records:
                out_fd.write(record.table_form2())

    def filter_by_function(self, function):

        filtered_records, filtered_out_records = self.filter_records_by_function(function)

        return (CollectionPrimer3(record_list=filtered_records, primer3_file=None, from_file=False),
                CollectionPrimer3(record_list=filtered_out_records, primer3_file=None, from_file=False))

    def filter_out_records_without_primers(self):

        def record_with_primers(record):
            if record.primer_pair_count > 0:
                return True
            return False

        return self.filter_by_function(record_with_primers)

    def remove_primers_with_gaps_in_pcr_product(self, min_gap_len=5):
        for record in self.records:
            record.remove_primers_with_gaps_in_pcr_product(min_gap_len=min_gap_len)

    def seq_ids(self):
        id_list = IdList()

        for record in self.records:
            id_list.append(record.id)

        return id_list

    def extract_records_by_ids(self, id_list):
        extracted_record_list = []
        for record in self.records:
            if record.id in id_list:
                extracted_record_list.append(record)

        return CollectionPrimer3(from_file=False, record_list=extracted_record_list)


