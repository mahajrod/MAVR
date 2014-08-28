#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class FastaAlternateReferenceMaker():
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html

    def correct_reference(self, gatk_dir, reference, new_reference, variants_vcf):

        os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --variant %s"
                  % (gatk_dir, reference, new_reference, variants_vcf))

    def restore_names(self, reference, new_reference, corrected_reference, sed_script="sed_script.scr", names_file="chr_name_lines.t"):
        os.system("grep '>' %s > %s" % (reference, names_file))
        with open(names_file, "r") as in_fd:
            with open(sed_script, "w") as out_fd:
                line_number = 1
                for line in in_fd:
                    out_fd.write("s/>%i$/%s/\n" % (line_number, line.strip()))
                    line_number += 1
        os.system("sed -f %s %s > %s" % (sed_script, new_reference, corrected_reference))

if __name__ == "__main__":
    workdir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/align_LAN210_reads_to_LAN210v0.9"
    os.chdir(workdir)
    FastaAlternateReferenceMaker = FastaAlternateReferenceMaker()
    FastaAlternateReferenceMaker.restore_names("../LAN210_v0.9m.fasta", "LAN210_v0.10m_bad_names.fasta", "LAN210_v0.10m.fasta")
