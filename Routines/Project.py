#!/usr/bin/env python
__author__ = 'mahajrod'

from collections import OrderedDict
from Routines.File import FileRoutines


class ProjectRoutines(FileRoutines):
    def __init__(self):
        FileRoutines.__init__(self)

    def initiate(self, workdir, project_name, species_list):
        """
        Structure of project folder:
            <workdir>/
                <project_name>/
                    <species1>
                        genome_denovo/
                            reads/
                                raw/
                                filtered/
                            fastqc/
                                raw/
                                filtered/
                            kmer/
                                raw/
                                filtered/
                            assemblies/
                            annotation/
                                protein_coding_genes/
                                repeats/
                                ncRNA/
                                    tRNA/
                                    rRNA/
                                    miRNA/
                                    other_ncRNA/
                            analysis/
                                orthologs/
                                families/
                                positive_selection/
                                specific_genes/
                        genome_ref_assisted/
                            reference/
                            reads/
                                raw/
                                filtered/
                            fastqc/
                                raw/
                                filtered/
                            kmer/
                                raw/
                                filtered/
                            SNPcall/
                            analysis/
                                population/
                        transcriptome/
                            reads/
                                raw/
                                filtered/
                            fastqc/
                                raw/
                                filtered/
                            kmer/
                                raw/
                                filtered/
                            assemblies/
                            annotation/
                            alignment/
                            analysis/
                            alignment/
                            expression/
                        analysis/

                    <species2>
                        ...
                    ...

                    <speciesN>
                        ...

        """
        species_dir_dict = {
                            "genome_denovo": {
                                "reads": {
                                    "raw": {},
                                    "filtered": {}, },
                                "fastqc": {
                                    "raw": {},
                                    "filtered": {}, },
                                "kmer": {
                                    "raw": {},
                                    "filtered": {}, },
                                "assemblies": {},
                                "annotation": {
                                    "protein_coding_genes": {},
                                    "repeats": {
                                        "repeatmasker": {},
                                        "TRF": {},
                                        "windowmasker": {}, },
                                    "ncRNA": {
                                        "tRNA": {},
                                        "rRNA": {},
                                        "miRNA": {},
                                        "other_ncRNA": {}, }, },
                                "analysis": {
                                    "orthologs": {},
                                    "families": {},
                                    "positive_selection": {},
                                    "specific_genes": {}, }
                                                },
                            "genome_ref_assisted": {
                                "reference": {},
                                "reads": {
                                    "raw": {},
                                    "filtered": {}, },
                                "fastqc": {
                                    "raw": {},
                                    "filtered": {}, },
                                "kmer": {
                                    "raw": {},
                                    "filtered": {}, },
                                "alignment": {},
                                "SNPcall": {
                                    "raw_vcf": {},
                                    "filtered_vcf": {},
                                    },
                                "analysis": {
                                    "population": {}, },
                                                    },
                            "transcriptome": {
                                "reads": {
                                    "raw": {},
                                    "filtered": {}, },
                                "fastqc": {
                                    "raw": {},
                                    "filtered": {}, },
                                "kmer": {
                                    "raw": {},
                                    "filtered": {}, },
                                "assemblies": {},
                                "annotation": {},
                                "alignment": {},
                                "analysis": {},
                                "expression": {},
                                              },
                            "analysis": {},
                            }

        manuscripts_dict = {
                            "plans": {},
                            "articles": {},
                            "drafts": {},
                             "figures": {},
                                            }

        project_dir_structure_dict = OrderedDict()
        project_dir_structure_dict[project_name] = OrderedDict()

        for species in [species_list] if isinstance(species_list, str) else species_list:
            project_dir_structure_dict[project_name][species] = species_dir_dict

        project_dir_structure_dict[project_name]["manuscripts"] = manuscripts_dict

        self.recursive_mkdir(project_dir_structure_dict, out_dir=workdir,
                             description_filename="DESCRIPTION", description_text=None,
                             readme_filename=None, readme_text=None)

