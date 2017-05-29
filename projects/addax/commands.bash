#!/usr/bin/env bash

#~/soft/MAVR/scripts/annotation/parallel_augustus.py -i K63g15d5e5R.scafSeq.gapclosed3.selected_repeat_classes.masked.fasta  -o addax_annotation -x ADDNAS -t 37 -s human -c ~/soft/augustus-3.2.1/co -p /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/db/pfam/Pfam-A.hmm -w /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/db/swissprot/swissprot --softmasking --hintsfile all.hints.gff --extrinsicCfgFile ~/soft/augustus-3.2.1/config/extrinsic/extrinsic.RM.EXNT.EXNS.cfg -a ~/soft/augustus-3.2.1/bin/

~/soft/MAVR/scripts/annotation/parallel_augustus.py -i K63g15d5e5R.scafSeq.gapclosed3.selected_repeat_classes.masked.fasta  \
                                                    -o addax_annotation \
                                                    -x ADDNAS \
                                                    -t 37 \
                                                    -s human \
                                                    -c ~/soft/augustus-3.2.1/config/ \
                                                    -p /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/db/pfam/Pfam-A.hmm \
                                                    -w /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/db/swissprot/swissprot \
                                                    --softmasking \
                                                    --hintsfile all.hints.gff \
                                                    --extrinsicCfgFile \
                                                    ~/soft/augustus-3.2.1/config/extrinsic/extrinsic.RM.EXNT.EXNS.cfg \
                                                    -a ~/soft/augustus-3.2.1/bin/