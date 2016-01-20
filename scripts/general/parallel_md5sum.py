


ls  */*.gz | xargs -P 6 -I NAME md5sum NAME > fastq.md5sum