

#dell
cd /mnt/guatemala/skliver/Boechera/Boechera_holboellii/WGA/masking/trf
for SP in `ls ../../genomes/`; do  mkdir ${SP}; cd ${SP}; ~/Soft/MAVR/scripts/repeat_masking/tandem_repeat_masking.py -i ../../../genomes/${SP}/${SP}.fasta -p ~/Soft/TRF/trf -o ${SP} -t 10 & cd ../;  done
