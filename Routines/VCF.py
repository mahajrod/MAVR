
from Routines import SequenceRoutines


class VCFRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)

    def combine_same_samples_vcfs(self, vcf_list, output, close_fd_after=False, extension_list=[".vcf", ]):

        output_fd = output if isinstance(output, file) else open(output, "w")

        vcf_files = sorted(self.make_list_of_path_to_files_by_extension(vcf_list, extension_list=extension_list,
                                                                        recursive=False, return_absolute_paths=True))

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
