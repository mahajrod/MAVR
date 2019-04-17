#!/usr/bin/env bash

RM_FILE="/home/mahajrod/Reference_genomes/homo_sapiens/repeatmasking/GCF_000001405.31_GRCh38.p5_rm.out"

for REPEAT_TYPE in DNA DNA? DNA/hAT DNA/hAT? DNA/hAT-Blackjack DNA/hAT-Charlie DNA/hAT-Tip100 DNA/hAT-Tip100? DNA/Merlin DNA/MuDR DNA/MULE-MuDR DNA/PiggyBac DNA/PiggyBac? DNA/TcMar DNA/TcMar? DNA/TcMar-Mariner DNA/TcMar-Pogo DNA/TcMar-Tc2 DNA/TcMar-Tigger LINE/CR1 LINE/Dong-R4 LINE/L1 LINE/L2 LINE/Penelope LINE/RTE-BovB LINE/RTE-X Low_complexity LTR LTR? LTR/ERV1 LTR/ERV1? LTR/ERVK LTR/ERVL LTR/ERVL? LTR/ERVL-MaLR LTR/Gypsy LTR/Gypsy? RC?/Helitron? RC/Helitron Retroposon RNA rRNA Satellite Satellite/acro Satellite/centr Satellite/telo scRNA Simple_repeat SINE? SINE/Alu SINE/Deu SINE/MIR SINE/tRNA snRNA srpRNA tRNA Unknown;
    do

    COUNT=`grep -c " ${REPEAT_TYPE} " /home/mahajrod/Reference_genomes/homo_sapiens/repeat_masking/GCF_000001405.31_GRCh38.p5_rm.out`
    echo -e "${REPEAT_TYPE}\t${COUNT}"
    done