#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
from collections import OrderedDict

from Bio import SeqIO


from Routines import SequenceRoutines, FileRoutines
from CustomCollections.GeneralCollections import IdList, IdSet


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use sys.stdout.write instead of print to safe write to stdout from multiple threads
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    os.system(exe_string)


class Tool():

    def __init__(self, cmd, path="", max_threads=4, jar_path=None, jar=None,
                 max_memory="500m", timelog="tool_time.log"):
        self.path = FileRoutines.check_path(path)
        self.cmd = cmd
        self.threads = max_threads
        #print(jar_path)
        self.jar_path = FileRoutines.check_path(jar_path) if jar_path else None
        self.jar = jar
        self.max_memory = max_memory
        self.timelog = timelog

    def execute(self, options="", cmd=None, capture_output=False):
        command = cmd if cmd is not None else self.cmd

        exe_string = (FileRoutines.check_path(self.path) if self.path else "") + command + " " + options

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        exe_string = "time -f 'Time\\t%%E real,\\t%%U user,\\t%%S sys' -a -o %s %s" % (self.timelog, exe_string) if self.timelog else exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None):
        command = cmd if cmd is not None else self.cmd
        exe_string_list = [(FileRoutines.check_path(self.path) if self.path else "") + command + " " + options
                           for options in options_list]

        with open("exe_list.t", "a") as exe_fd:
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
        FileRoutines.save_mkdir(output_dir)
        out_prefix = FileRoutines.split_filename(input_fasta)[1] if output_prefix is None else output_prefix
        sequence_dict = SeqIO.index_db("temp.idx", input_fasta, "fasta")

        split_index = 1
        records_written = 0
        record_ids_list = list(sequence_dict.keys())
        number_of_records = len(record_ids_list)

        num_of_recs = int(number_of_records/num_of_files) + 1 if num_of_files else num_of_recs_per_file
        while (records_written + num_of_recs) <= number_of_records:

            SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict,
                                                                record_ids_list[records_written:records_written+num_of_recs]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
            split_index += 1
            records_written += num_of_recs

        if records_written != number_of_records:
            SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict,
                                                                record_ids_list[records_written:]),
                        "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

        os.remove("temp.idx")

    def split_fasta_by_seq_len(self, input_fasta, output_dir, max_len_per_file=None, output_prefix=None):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        FileRoutines.save_mkdir(output_dir)

        out_prefix = FileRoutines.split_filename(input_fasta)[1] if output_prefix is None else output_prefix
        sequence_dict = SeqIO.index_db("temp.idx", input_fasta, "fasta")
        length = 0

        for record_id in sequence_dict:
            length += len(sequence_dict[record_id].seq)
            
        max_len = max_len_per_file if max_len_per_file else int(length / self.threads)

        split_index = 1
        id_list = []
        total_length = 0

        for record_id in sequence_dict:
            record_length = len(sequence_dict[record_id].seq)
            if record_length >= max_len:
                SeqIO.write(sequence_dict[record_id],
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

            elif total_length + record_length > max_len:
                SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = record_length
                id_list = [record_id]

            elif total_length + record_length == max_len:
                id_list.append(record_id)
                SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")
                total_length = 0
                id_list = []

            elif total_length + record_length < max_len:
                id_list.append(record_id)
                total_length += record_length
                continue

            split_index += 1

        if id_list:
            SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict, id_list),
                            "%s/%s_%i.fasta" % (output_dir, out_prefix, split_index), format="fasta")

        os.remove("temp.idx")

    @staticmethod
    def extract_common_sequences(list_of_files_with_sequences_of_samples, list_of_names_of_samples,
                                 output_dir, separator="_", format="fasta"):

        def generator_with_id_correction(samples_seq_dict, common_record_id):
            for sample in samples_seq_dict:
                record = samples_seq_dict[sample][common_record_id]
                record.id = "%s%s%s" % (sample, separator, record.id)
                yield record

        FileRoutines.save_mkdir(output_dir)
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
                        "%s%s.%s" % (FileRoutines.check_path(output_dir), common_id, format),
                        format=format)
        for i in range(0, index):
            os.remove("tmp_%i.idx" % i)

    @staticmethod
    def extract_ids_from_file(input_file, output_file=None, header=False, column_separator="\t",
                              comments_prefix="#", column_number=None):
        id_list = IdList()
        id_list.read(input_file, column_separator=column_separator, comments_prefix=comments_prefix,
                     column_number=column_number, header=header)
        if output_file:
            id_list.write(output_file, header=header)
        return id_list

    @staticmethod
    def intersect_ids_from_files(files_with_ids_from_white_list, files_with_ids_from_black_list, result_file=None):
        white_set = IdSet()
        black_set = IdSet()

        for filename in files_with_ids_from_white_list:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#")
            white_set = white_set | id_set

        for filename in files_with_ids_from_black_list:
            id_set = IdSet()
            id_set.read(filename, comments_prefix="#")
            black_set = black_set | id_set

        final_set = IdSet(white_set - black_set)

        final_set.write(result_file if result_file else sys.stdout)

    @staticmethod
    def intersect_ids(list_of_white_id_lists, list_of_black_id_list):
        white_set = IdSet()
        black_set = IdSet()

        for id_list in list_of_white_id_lists:
            white_set = white_set | IdSet(id_list)

        for id_list in list_of_black_id_list:
            black_set = black_set | IdSet(id_list)
        final_set = IdSet(white_set - black_set)

        return final_set


class JavaTool(Tool):

    def __init__(self, jar, java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):

        Tool.__init__(self, "java", path=java_path, max_threads=max_threads,
                      jar_path=jar_path, jar=jar, max_memory=max_memory, timelog=timelog)

    def execute(self, options="", cmd=None, capture_output=False):
        command = cmd if cmd is not None else ""

        java_string = "java"
        java_string += " -Xmx%s" if self.max_memory else ""
        #print (self.jar_path)
        java_string += " -jar %s%s" % (self.jar_path, self.jar)
        java_string += " %s" % command
        java_string += " %s" % options

        exe_string = (FileRoutines.check_path(self.path) if self.path else "") + java_string

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        exe_string = "time -f 'Time\\t%%E real,\\t%%U user,\\t%%S sys' -a -o %s %s" % (self.timelog, exe_string) if self.timelog else exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

