#--------------------------------Variables---------------------------------------

WORKDIR="/mnt/peru/skliver/Boechera/Boechera_holboellii/genome/annotations/final/"

REPEATMASKING_DIR="${WORKDIR}/repeat_masking/"

EXONERATE_DIR="${WORKDIR}/exonerate/"
EXONERATE_PEP_DIR="${EXONERATE_DIR}/protein/"
EXONERATE_TRANSCRIPT_DIR="${EXONERATE_DIR}/transcript"


AUGUSTUS_DIR="${WORKDIR}/augustus/"



#--------------------------------------------------------------------------------


#-------------------------------------Script-------------------------------------

#repeatmasking

cd ${REPEATMASKING_DIR}

~/Soft/MAVR/scripts/repeat_masking/convert_rm_out_to_gff.py -i gap_closed_assembly.250.fasta.out -p gap_closed_assembly.250.fasta.converted

bedtools maskfasta -fi gap_closed_assembly.250.fasta -bed gap_closed_assembly.250.selected_repeat_classes.gff  -soft -fo gap_closed_assembly.250.softmasked.fasta


# prepare hints

cd ${EXONERATE_PEP_DIR}

for DIR in `ls`; do echo ${DIR}; cd ${DIR}; ~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i ./out -o ${DIR}&  cd ../;  done

for DIR in `ls`; do echo ${DIR}; cd ${DIR}; ~/Soft/MAVR/scripts/annotation/prepare_hints_from_exonerate_target_output.py -i ${DIR}.target.gff -o ${DIR}.hints -e ~/Soft/augustus-3.2.1/scripts/ &  cd ../;  done


cd ${EXONERATE_TRANSCRIPT_DIR}

for DIR in `ls`; do echo ${DIR}; cd ${DIR}; ~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i ./out -o ${DIR}&  cd ../;  done

for DIR in `ls`; do echo ${DIR}; cd ${DIR}; ~/Soft/MAVR/scripts/annotation/prepare_hints_from_exonerate_target_output.py -i ${DIR}.target.gff -o ${DIR}.hints -e ~/Soft/MAVR/scripts/annotation/augustus/ --top_hits_CDS_part_cutoff 0 --secondary_hits_CDS_part_cutoff 3 --source_for_top_hits EXNTRT --source_for_secondary EXNTRS --top_hits_priority 40 --secondary_hits_priority 10 &  cd ../;  done


cd ${AUGUSTUS_DIR}

cp ${EXONERATE_PEP_DIR}/*/*hints.gff ${EXONERATE_TRANSCRIPT_DIR}/*/*hints.gff ./

cat *.hints.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > all.hints.gff

ln /mnt/peru/skliver/Boechera/Boechera_holboellii/genome/annotations/final/repeat_masking/gap_closed_assembly.250.softmasked.fasta ./


~/Soft/MAVR/scripts/annotation/parallel_augustus.py -i gap_closed_assembly.250.softmasked.fasta -o boechera_retrofracta -t 40 -x BOERET -s arabidopsis -c ~/Soft/augustus-3.2.1/config/ -p /mnt/guatemala/skliver/data/Pfam/Pfam-A.hmm -w /mnt/guatemala/skliver/data/SwissProt/swissprot --softmasking --hintsfile all.hints.gff  --extrinsicCfgFile ~/Soft/augustus-3.2.1/config/extrinsic/extrinsic.RM.EXNT.EXNS.RNASEQ.EXNTTRT.EXNTTRS.cfg -a ~/Soft/augustus-3.2.1/bin/
#--------------------------------------------------------------------------------

screen