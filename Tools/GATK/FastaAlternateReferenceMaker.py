#!/usr/bin/env python
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from Tools.Abstract import JavaTool

from Routines import AnnotationsRoutines
from CustomCollections.GeneralCollections import SynDict


class FastaAlternateReferenceMaker(JavaTool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html

    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    @staticmethod
    def parse_options(reference, new_reference, variants_vcf, raw_seq_per_line=False, vcf_with_masking=None,
                      override_vcf_by_mask=None, use_ambiguous_nuccleotides=None, interval_list=None):
        options = " -R %s" % reference
        options += " -o %s" % new_reference
        options += " --variant %s" % variants_vcf
        options += " --rawOnelineSeq" if raw_seq_per_line else ""
        options += " --snpmask %s" % vcf_with_masking if vcf_with_masking else ""
        options += " --snpmaskPriority" if override_vcf_by_mask else ""
        options += " --use_IUPAC" if use_ambiguous_nuccleotides else ""

        if interval_list:
            if isinstance(interval_list, str):
                options += " -L %s" % interval_list
            else:
                for entry in interval_list:
                    options += " -L %s" % entry

        return options

    def correct_reference(self,
                          reference,
                          new_reference,
                          variants_vcf,
                          raw_seq_per_line=False,
                          vcf_with_masking=None,
                          override_vcf_by_mask=None,
                          use_ambiguous_nuccleotides=None,
                          interval_list=None):

        options = self.parse_options(reference,
                                     new_reference,
                                     variants_vcf,
                                     raw_seq_per_line=raw_seq_per_line,
                                     vcf_with_masking=vcf_with_masking,
                                     override_vcf_by_mask=override_vcf_by_mask,
                                     use_ambiguous_nuccleotides=use_ambiguous_nuccleotides,
                                     interval_list=interval_list)

        self.execute(options=options)

        #os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T FastaAlternateReferenceMaker -o %s --variant %s"
        #          % (gatk_dir, reference, new_reference, variants_vcf))

    def correct_regions_from_gff(self,
                                 reference,
                                 variants_vcf,
                                 gff_file,
                                 output_prefix=None,
                                 feature_type_list=["CDS"],
                                 unification_key="Parent",
                                 #raw_seq_per_line=False,
                                 vcf_with_masking=None,
                                 override_vcf_by_mask=None,
                                 use_ambiguous_nuccleotides=None):

        feature_dict = AnnotationsRoutines.get_feature_dict(gff_file,
                                                            output_prefix=output_prefix,
                                                            feature_type_list=feature_type_list,
                                                            unification_key=unification_key)
        region_file = "%s.coordinates_only.list" % output_prefix

        raw_regions = "%s.raw.seq" % output_prefix
        final_regions = "%s.fasta" % output_prefix

        self.correct_reference(reference,
                               raw_regions,
                               variants_vcf,
                               raw_seq_per_line=True,
                               vcf_with_masking=vcf_with_masking,
                               override_vcf_by_mask=override_vcf_by_mask,
                               use_ambiguous_nuccleotides=use_ambiguous_nuccleotides,
                               interval_list=region_file)

        region_with_frameshift = SynDict()

        def new_regions_generator():
            with open(raw_regions, "r") as in_fd:
                for region_id in feature_dict:
                    seq = ""
                    for i in range(0, len(feature_dict[region_id])):
                        seq_fragment = in_fd.readline().strip()
                        if ((int(feature_dict[region_id][i][2]) - int(feature_dict[region_id][i][1]) + 1) - len(seq_fragment)) % 3 != 0:
                            if region_id not in region_with_frameshift:
                                region_with_frameshift[region_id] = [i]
                            else:
                                region_with_frameshift[region_id].append(i)
                        seq += seq_fragment
                    yield SeqRecord(seq=Seq(seq) if feature_dict[region_id][0][3] == "+" else Seq(seq).reverse_complement(),
                                    id=region_id,
                                    description="")

        SeqIO.write(new_regions_generator(), final_regions, format="fasta")

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
