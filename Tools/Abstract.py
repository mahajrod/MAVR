#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
from collections import OrderedDict

from Bio import SeqIO

from Routines.File import check_path, split_filename, save_mkdir
from Routines.Sequence import record_by_id_generator


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use sys.stdout.write instead of print to safe write to stdout from multiple threads
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    os.system(exe_string)


class Tool():

    def __init__(self, cmd, path="", max_threads=4, jar_path=None, jar=None, max_memory="500m"):
        self.path = check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        self.jar_path = check_path(jar_path) if jar_path else None
        self.jar = jar
        self.max_memory = max_memory

    def execute(self, options="", cmd=None, capture_output=False):

        command = cmd if cmd is not None else self.cmd
        exe_string = check_path(self.path) + command + " " + options
        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None):
        command = cmd if cmd is not None else self.cmd
        exe_string_list = [check_path(self.path) + command + " " + options for options in options_list]
        with open("exe_list.t", "w") as exe_fd:
            for entry in exe_string_list:
                exe_fd.write("%s\n" % entry)
        process_pool = mp.Pool(threads if threads else self.threads)

        results = process_pool.map(execute, exe_string_list)
        #process_pool.close()
        return results if capture_output else None

    @staticmethod
    def split_fasta(input_fasta, output_dir, num_of_recs_per_file=None, num_of_files=None, output_prefix=None):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        save_mkdir(output_dir)
        out_prefix = split_filename(input_fasta)[1] if output_prefix is None else output_prefix
        sequence_dict = SeqIO.index_db("temp.idx", input_fasta, "fasta")

        split_index = 1
        records_written = 0
        record_ids_list = list(sequence_dict.keys())
        number_of_records = len(record_ids_list)

        num_of_recs = int(number_of_records/num_of_files) + 1 if num_of_files else num_of_recs_per_file
        while (records_written + num_of_recs) <= number_of_records:

            SeqIO.write(record_by_id_generator(sequence_dict,
                                               record_ids_list[records_written:records_written+num_of_recs]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
            split_index += 1
            records_written += num_of_recs

        if records_written != number_of_records:
            SeqIO.write(record_by_id_generator(sequence_dict,
                                               record_ids_list[records_written:]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
        os.remove("temp.idx")

    @staticmethod
    def split_fasta_by_seq_len(input_fasta, output_dir, max_len_per_file=100000, output_prefix=None):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        save_mkdir(output_dir)

        out_prefix = split_filename(input_fasta)[1] if output_prefix is None else output_prefix
        sequence_dict = SeqIO.index_db("temp.idx", input_fasta, "fasta")

        split_index = 1
        id_list = []
        total_length = 0

        for record_id in sequence_dict:
            record_length = len(sequence_dict[record_id].seq)
            if record_length >= max_len_per_file:
                SeqIO.write(sequence_dict[record_id],
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

            elif total_length + record_length > max_len_per_file:
                SeqIO.write(record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = record_length
                id_list = [record_id]

            elif total_length + record_length == max_len_per_file:
                id_list.append(record_id)
                SeqIO.write(record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = 0
                id_list = []

            elif total_length + record_length < max_len_per_file:
                id_list.append(record_id)
                total_length += record_length
                continue

            split_index += 1

        os.remove("temp.idx")

    @staticmethod
    def extract_common_sequences(list_of_files_with_sequences_of_samples, list_of_names_of_samples,
                                 output_dir, separator="_", format="fasta"):

        def generator_with_id_correction(samples_seq_dict, common_record_id):
            for sample in samples_seq_dict:
                record = samples_seq_dict[sample][common_record_id]
                record.id = "%s%s%s" % (sample, separator, record.id)
                yield record

        save_mkdir(output_dir)
        index = 0
        samples_seq_dict = OrderedDict()
        for filename, sample_name in zip(list_of_files_with_sequences_of_samples, list_of_names_of_samples):
            samples_seq_dict[sample_name] = SeqIO.index_db("tmp_%i.idx" % index, filename, format=format)
            index += 1
        common_sequence_ids = set(samples_seq_dict[list_of_names_of_samples[0]].keys())

        for sample_name in list_of_names_of_samples[1:]:
            common_sequence_ids = common_sequence_ids & set(samples_seq_dict[sample_name].keys())

        for common_id in common_sequence_ids:
            SeqIO.write(generator_with_id_correction(samples_seq_dict, common_id),
                        "%s%s.%s" % (check_path(output_dir), common_id, format),
                        format=format)
        for i in range(0, index):
            os.remove("tmp_%i.idx" % i)


class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory="1g"):

        Tool.__init__(self, "java", path=check_path(java_path), max_threads=max_threads,
                      jar_path=check_path(jar_path),
                      jar=jar, max_memory=max_memory)

    def execute(self, options="", cmd=None):
        command = cmd if cmd is not None else ""
        exe_string = check_path(self.path) + "java -Xmx%s -jar %s%s %s %s" % (self.max_memory,
                                                                              check_path(self.jar_path),
                                                                              self.jar, command, options)
        sys.stdout.write("Executing:\n\t%s\n" % exe_string)

        os.system(exe_string)


if __name__ == "__main__":
    os.chdir("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/")
    Tool.extract_common_sequences(["LAN210_v0.10m_selected_proteins.fasta", "S288C_R64_selected_proteins.fasta"],
                                  ["LAN210", "S288C"], "proteins")
