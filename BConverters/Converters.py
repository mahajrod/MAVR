__author__ = 'mahajrod'

import re
from Bio import Phylo


def convert_tree(input_file, input_filetype, output_file, output_filetype):
    #tree = Phylo.read(input_file, input_filetype)
    Phylo.convert(input_file, input_filetype, output_file, output_filetype)
    #Phylo.write(tree, output_file, output_filetype)


def sam2gtf(input_file, output_file, start_position=1, source="custom", scaffold=19):
    #TODO: works only for transcript alignments!!!!!!!!!!!!
    isoform = ["a", "b", "c", "d"]
    shift = start_position - 1
    isoform_index = 0
    with open(input_file, "r") as in_fd:
        with open(output_file, "w") as out_fd:
            for line in in_fd:
                if line[0] == "@":
                    continue
                sam_line = line.strip().split("\t")
                id = sam_line[0]
                region = sam_line[2]
                start = int(sam_line[3])
                quality = sam_line[4]
                cigar_string = sam_line[5]
                splited = re.split("M|N", cigar_string)[:-1]
                print(splited)
                segments_lengthes = list(map(lambda x: int(x), splited))
                isoform_name = "nxf1-" + isoform[isoform_index]
                isoform_index += 1
                transcript_line = '19\t%s\ttranscript\t%i\t%i\t.\t+\t.\tgene_id "NXF1"; gene_version "1"; transcript_id "%s"; transcript_version "1"; gene_name "nxf-1"; gene_source "custom"; gene_biotype "protein_coding"; transcript_name "%s";\n' % (source, start + shift, shift + sum(segments_lengthes), isoform_name, isoform_name)
                out_fd.write(transcript_line)
                for i in range(0, len(segments_lengthes)):
                    exon_start = shift + start + sum(segments_lengthes[:i])
                    if i % 2 == 0:
                        exon_line = '19\t%s\texon\t%i\t%i\t.\t+\t.\tgene_id "NXF1"; gene_version "1"; transcript_id "%s"; transcript_version "1"; exon_number "%i"; gene_name "nxf-1"; gene_source "custom"; gene_biotype "protein_coding"; transcript_name "%s";\n' % (source, exon_start, exon_start + segments_lengthes[i] - 1, isoform_name, i / 2, isoform_name)
                        out_fd.write(exon_line)



    #skip_header