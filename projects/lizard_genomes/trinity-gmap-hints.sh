~/Soft/MAVR/scripts/transcriptome/alignment/filter_gmap_alignments.py -i Trinita.gff -o Trinity.filtered.coverage90.identity90.gff -c 90 -v -d 90
~/Soft/MAVR/scripts/annotation/augustus/gff2hints.pl --priority=120  --CDSpart_cutoff=0 --source=RNASEQ --in=Trinity.transcript.gff --out=Trinity.hints.gff --hintid=CDS --transcriptfeaturetype=mRNA
~/Soft/MAVR/scripts/annotation/augustus/gff2hints.pl --priority=120  --CDSpart_cutoff=0 --source=RNASEQ --in=Trinity.filtered.coverage90.identity90.gff --out=Trinity.hints.gff --hintid=CDS --transcriptfeaturetype=mRNA
