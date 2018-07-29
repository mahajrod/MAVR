__author__ = 'mahajrod'
import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Tools.Abstract import Tool

#TODO: finish code
def examine_gff(gff_file):
    examiner = GFFExaminer()
    in_handle = open(gff_file)
    pprint.pprint(examiner.available_limits(in_handle))
    print("")
    in_handle.close()


def parse_gff(gff_file, limit_info={}):
    examine_gff(gff_file)
    file_handler = open(gff_file)
    if limit_info:
        gff_ann = list(GFF.parse(file_handler, limit_info=limit_info))
    else:
        gff_ann = list(GFF.parse(file_handler))
    file_handler.close()
    #print gff_ann
    return gff_ann


def make_glimmer_training_data(seq_file, seq_index, gff_file, out_transcript_file,
                               out_exon_coordinates, seq_filetype="fasta"):
    fasta_dict = SeqIO.index_db(seq_index, [seq_file], seq_filetype)

    def generator(gff_file, fasta_dict, exon_fd):
        with open(gff_file, "r") as gff_fd:
            for record in GFF.parse(gff_fd, target_lines=100000):
                for feature in record.features:
                    #print (feature.type)
                    if feature.type == "transcript":
                        print(feature)
                        print(feature.sub_features)
                        exon_fd.write(feature.id + "\n")
                        exon_fd.write(str(feature.location))
                        exon_fd.write("\n")
                        exon_fd.write(str(feature.sub_features) + "\n")

                    if feature.type == "gene":
                        feature_record = feature.extract(fasta_dict[record.id])
                        feature_record.id = feature.qualifiers["gene_id"][0]
                        feature_record.description = ""
                        yield feature_record
    exon_fd = open(out_exon_coordinates, "w")
    SeqIO.write(generator(gff_file, fasta_dict, exon_fd), out_transcript_file, "fasta")
    exon_fd.close()


if __name__ == "__main__":

    make_glimmer_training_data("/run/media/mahajrod/Data/ensemble_75/fasta/takifugu_rubripes/dna/Takifugu_rubripes.FUGU4.75.dna.toplevel.fa",
                               "/run/media/mahajrod/Data/ensemble_75/fasta/takifugu_rubripes/dna/Takifugu_rubripes.FUGU4.75.dna.toplevel.index",
                               "/run/media/mahajrod/Data/ensemble_75/gtf/takifugu_rubripes/Takifugu_rubripes.FUGU4.75.gtf",
                               "/run/media/mahajrod/Data/ensemble_75/fasta/takifugu_rubripes/dna/Takifugu_rubripes.FUGU4.75.transcripts.fa",
                               "/run/media/mahajrod/Data/ensemble_75/fasta/takifugu_rubripes/dna/Takifugu_rubripes.FUGU4.75.exons.coordinates",
                               seq_filetype="fasta")
