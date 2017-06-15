__author__ = 'mahajrod'
from Tools.Abstract import Tool


class Cutadapt(Tool):

    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "cutadapt", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(forward_file, output_prefix, reverse_file=None, format="fastq",
                      three_prime_adapter_list=None, five_prime_adapter_list=None, anyway_adapter_list=None,
                      max_number_of_adapters_per_read=None, trim_Ns_on_read_end=False, min_read_length_after_trimming=None):

        three_prime_adapters = [three_prime_adapter_list] if isinstance(three_prime_adapter_list, str) else three_prime_adapter_list
        five_prime_adapters = [five_prime_adapter_list] if isinstance(five_prime_adapter_list, str) else five_prime_adapter_list
        anyway_adapters = [anyway_adapter_list] if isinstance(anyway_adapter_list, str) else anyway_adapter_list

        options = ""
        for three_adapter in three_prime_adapters:
            options = " -a %s" % three_adapter

        for five_adapter in five_prime_adapters:
            options = " -g %s" % five_adapter

        for anyway_adapter in anyway_adapters:
            options = " -b %s" % anyway_adapter

        options += " -n %i" % max_number_of_adapters_per_read if max_number_of_adapters_per_read else ""
        options += " --trim-n" if trim_Ns_on_read_end else ""
        options += " -m %i" % min_read_length_after_trimming if min_read_length_after_trimming else ""

        if reverse_file:
            filtered_forward = "%s_1.%s" % (output_prefix, format)
            filtered_reverse = "%s_2.%s" % (output_prefix, format)

            options += " -o %s" % filtered_forward
            options += " -p %s" % filtered_reverse
            options += " %s %s" % (forward_file, reverse_file)
        else:
            filtered_se = "%s.se.%s" % (output_prefix, format)
            options += " -o %s" % filtered_se
            options += " %s" % forward_file

        return options

    def filter(self, forward_file, output_prefix, reverse_file=None, format="fastq",
               three_prime_adapter_list=None, five_prime_adapter_list=None, anyway_adapter_list=None,
               max_number_of_adapters_per_read=None, trim_Ns_on_read_end=False, min_read_length_after_trimming=None):

        options = self.parse_options(forward_file, output_prefix, reverse_file=reverse_file,
                                     format=format,
                                     three_prime_adapter_list=three_prime_adapter_list,
                                     five_prime_adapter_list=five_prime_adapter_list,
                                     anyway_adapter_list=anyway_adapter_list,
                                     max_number_of_adapters_per_read=max_number_of_adapters_per_read,
                                     trim_Ns_on_read_end=trim_Ns_on_read_end,
                                     min_read_length_after_trimming=min_read_length_after_trimming)

        self.execute(options=options)








