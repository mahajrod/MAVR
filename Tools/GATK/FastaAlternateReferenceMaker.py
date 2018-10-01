#!/usr/bin/env python
import os
from Tools.Abstract import JavaTool


class FastaAlternateReferenceMaker(JavaTool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html

    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def correct_reference(self, reference, new_reference, variants_vcf):

        options = " -R %s" % reference
        options += " -o %s" % new_reference
        options += " --variant %s" % variants_vcf

        self.execute(options=options)

        #os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --variant %s"
        #          % (gatk_dir, reference, new_reference, variants_vcf))

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
