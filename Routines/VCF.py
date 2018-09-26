import os

from collections import OrderedDict

from CustomCollections.GeneralCollections import IdList, SynDict
from Routines.Sequence import SequenceRoutines


class VCFRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)

    def combine_same_samples_vcfs(self, vcf_list, output, close_fd_after=False, extension_list=[".vcf", ],
                                  order_vcf_files=False, sort=False):

        output_fd = output if isinstance(output, file) else open(output, "w")

        vcf_files = self.make_list_of_path_to_files_by_extension(vcf_list, extension_list=extension_list,
                                                                 recursive=False, return_absolute_paths=True)
        if order_vcf_files:
            vcf_files.sort()

        if sort:
            unsorted_file = "%s.unsorted.tmp" % output

            #string = "cat %s > %s" % (vcf_files[0], unsorted_file)
            #os.system(string)
            #with open(unsorted_file, "w") as out_fd:
            #    pass

            for filename in vcf_files[1:]:
                string = "sed -n '/^[^#]/p' %s >> %s" % (filename, unsorted_file)
                print(string)
                os.system(string)

            sorting_string = "(sed '/^[^#]/Q' %s; sort -k1,1 -k2,2n %s) > %s" % (vcf_files[0],
                                                                                 unsorted_file,
                                                                                 output)
            print(sorting_string)

            os.system(sorting_string)

        else:
            with open(vcf_files[0], "r") as in_fd:
                for line in in_fd:
                    output_fd.write(line)

            if len(vcf_files) > 1:
                for filename in vcf_files[1:]:
                    with open(filename, "r") as in_fd:
                        for line in in_fd:
                            if line[0] == "#":
                                continue
                            else:
                                output_fd.write(line)
                                break
                        for line in in_fd:
                            output_fd.write(line)

        if not isinstance(output, file) or close_fd_after:
            output_fd.close()

    def check_gvcf_integrity(self, gvcf_file, output_prefix, reference=None, length_dict=None, parsing_mode="parse"):
        len_dict = length_dict if length_dict else self.get_lengths(record_dict=self.parse_seq_file(reference,
                                                                                                    mode=parsing_mode),
                                                                    out_file=None,
                                                                    close_after_if_file_object=False)

        scaffold_dict = OrderedDict()

        with self.metaopen(gvcf_file, "r") as gvcf_fd:
            prev_scaffold = ""

            for line in gvcf_fd:
                #print line
                if line[0] == "#":
                    continue

                line_list = line.split("\t")
                scaffold = line_list[0]
                start = int(line_list[1])
                format = line_list[7].split(";")

                if (len(format) == 1) and (format[0][0:3] == "END"):
                    end = int(format[0].split("=")[1])
                else:
                    end = start + len(line_list[3]) - 1
                #print line_list
                #print scaffold, start, end, format

                if scaffold not in scaffold_dict:
                    scaffold_dict[scaffold] = []

                if scaffold != prev_scaffold:
                    scaffold_dict[scaffold].append([start, end])
                else:
                    #print scaffold_dict[scaffold][-1][1]
                    if scaffold_dict[scaffold][-1][1] + 1 >= start:
                        scaffold_dict[scaffold][-1][1] = max(end, scaffold_dict[scaffold][-1][1])
                    else:
                        scaffold_dict[scaffold].append([start, end])
                prev_scaffold = scaffold

        complete_scaffolds = IdList()
        fragmented_scaffolds = IdList()
        scaffolds_with_absent_fragments = IdList()

        with open("%s.scaffold_regions" % output_prefix, "w") as scaf_reg_fd:

            for scaffold in scaffold_dict:
                if len(scaffold_dict[scaffold]) > 1:
                    fragmented_scaffolds.append(scaffold)

                scaffold_length = sum(map(lambda s: s[1] - s[0] + 1, scaffold_dict[scaffold]))
                if scaffold_length != len_dict[scaffold]:
                    scaffolds_with_absent_fragments.append(scaffold)
                else:
                    complete_scaffolds.append(scaffold)
                scaf_reg_fd.write("%s\t%s\n" % (scaffold, ",".join(map(lambda s: "-".join(map(str,s)), scaffold_dict[scaffold]))))

        complete_scaffolds.write("%s.complete_scaffolds" % output_prefix)
        fragmented_scaffolds.write("%s.fragmented_scaffolds" % output_prefix)
        scaffolds_with_absent_fragments.write("%s.scaffolds_with_absent_fragments" % output_prefix)