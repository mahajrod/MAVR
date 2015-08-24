#!/usr/bin/env bash
#!/usr/bin/env bash
WORKDIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/"

GENE_ALIAS_FILE=${WORKDIR}"genes_rev.syn"
INPUT_VCF_DIR=${WORKDIR}"annotated_vcf/"
#VCF_WITH_SELECTED_SNPS_DIR=${WORKDIR}"vcf_with_selected_snps/"
SNPEFF_ALL_INFO_DIR=${WORKDIR}"snpeff_all_info/"
SNPEFF_SELECTED_SNPS_INFO_DIR=${WORKDIR}"snpeff_selected_snps_info/"
SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR=${WORKDIR}"snpeff_selected_snps_not_silent_info/"
DESIRED_GENES_ID_FILE=${WORKDIR}genes.ids
SUMMARY_FILE=${WORKDIR}"snpeff_selected_snps_not_silent_info_summary.info"
SUMMARY_FILE_COMMON_GENE_NAME=${WORKDIR}"snpeff_selected_snps_not_silent_info_summary_common_gene_name.info"

SUMMARY_FILE_SINGLE_LETTER_AA=${WORKDIR}"snpeff_selected_snps_not_silent_info_summary_single_letter_aa.info"
SUMMARY_FILE_COMMON_GENE_NAME_SINGLE_LETTER_AA=${WORKDIR}"snpeff_selected_snps_not_silent_info_summary_common_gene_name_single_letter_aa.info"

GET_SNPEFF_INFO_SCRIPT="/home/mahajrod/Genetics/MACE/scripts/extract_snp_eff_info.py"
SUMMARY_SCRIPT=~/Genetics/MAVR/scripts/annotation/get_summary_table.py

cd ${WORKDIR}
#mkdir -p ${VCF_WITH_SELECTED_SNPS_DIR} ${SNPEFF_ALL_INFO_DIR} ${SNPEFF_SELECTED_SNPS_INFO_DIR} ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR}
mkdir -p ${SNPEFF_ALL_INFO_DIR} ${SNPEFF_SELECTED_SNPS_INFO_DIR} ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR}

for SAMPLE in `ls ${INPUT_VCF_DIR}/*.vcf`;
    do

    echo 'Handling ${SAMPLE}'
    SAMPLE_NAME=`basename ${SAMPLE} .vcf`

    SNPEFF_ALL_INFO_FILE=${SNPEFF_ALL_INFO_DIR}${SAMPLE_NAME}'_all.snpeff'
    SNPEFF_SELECTED_SNPS_INFO_FILE=${SNPEFF_SELECTED_SNPS_INFO_DIR}${SAMPLE_NAME}'_selected_snps.snpeff'
    SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_FILE=${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR}${SAMPLE_NAME}'_selected_snps_not_silent.snpeff'

    ${GET_SNPEFF_INFO_SCRIPT} -i ${SAMPLE} -o ${SNPEFF_ALL_INFO_FILE}

    head -n 1 ${SNPEFF_ALL_INFO_FILE} > ${SNPEFF_SELECTED_SNPS_INFO_FILE}
    grep -f ${DESIRED_GENES_ID_FILE} ${SNPEFF_ALL_INFO_FILE} >> ${SNPEFF_SELECTED_SNPS_INFO_FILE}

    grep -vP '\tSILENT\t' ${SNPEFF_SELECTED_SNPS_INFO_FILE} > ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_FILE}
    done

${SUMMARY_SCRIPT} -i ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR} -s "_annotated_selected_snps_not_silent" \
                  -o ${SUMMARY_FILE}
${SUMMARY_SCRIPT} -i ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR} -s "_annotated_selected_snps_not_silent" \
                  -o ${SUMMARY_FILE_COMMON_GENE_NAME} -g ${GENE_ALIAS_FILE}

${SUMMARY_SCRIPT} -i ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR} -s "_annotated_selected_snps_not_silent" \
                  -o ${SUMMARY_FILE_SINGLE_LETTER_AA} -r -c
${SUMMARY_SCRIPT} -i ${SNPEFF_SELECTED_SNPS_INFO_NOT_SILENT_DIR} -s "_annotated_selected_snps_not_silent" \
                  -o ${SUMMARY_FILE_COMMON_GENE_NAME_SINGLE_LETTER_AA} -g ${GENE_ALIAS_FILE} -r -c