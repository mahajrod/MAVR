#!/usr/bin/env python
import os
import multiprocessing as mp
from subprocess import PIPE, Popen

from Bio import SeqIO

from Routines.File import check_path, split_filename
from Routines.Sequence import record_by_id_generator


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
        print("Executing:\n\t%s" % exe_string)
        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            os.system(exe_string)
            return None

    def parallel_execute(self, options_list, cmd=None, capture_output=False, threads=None):
        command = cmd if cmd is not None else self.cmd
        exe_string_list = [check_path(self.path) + command + " " + options for options in options_list]

        process_pool = mp.Pool(threads if threads else self.threads)

        def execute(exe_string):
            print("Executing:\n\t%s" % exe_string)
            os.system(exe_string)

        results = process_pool.map(execute, exe_string_list)
        process_pool.close()
        return results if capture_output else None

    @staticmethod
    def split_fasta(input_fasta, output_dir, num_of_recs_per_file, num_of_files=None, output_prefix=None):
        """
        by default splits input files into files with num_of_recs_per_file.
        if num_of_files is set num_of_recs_per_file is ignored.
        """
        try:
            os.mkdir(output_dir)
        except OSError:
            pass
        out_prefix = split_filename(input_fasta)[1] if output_prefix is not None else output_prefix
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
        print("Executing:\n\t%s" % exe_string)

        os.system(exe_string)

