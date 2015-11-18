#!/usr/bin/env bash

#cd ~/workdir/cats/caracal
#~/soft/augustus-3.2.1/bin/augustus --species=Felis_catus --AUGUSTUS_CONFIG_PATH=/home/skliver/soft/augustus-3.2.1/config/ \
#                                   --strand=both --genemodel=complete \
#                                   --gff3=on /home/skliver/data/genomes/caracal/final.assembly.fasta > caracal_augustus_predictions.gff
#
#cd ~/workdir/cats/fishing_cat
#~/soft/augustus-3.2.1/bin/augustus --species=Felis_catus --AUGUSTUS_CONFIG_PATH=/home/skliver/soft/augustus-3.2.1/config/ \
#                                   --strand=both --genemodel=complete \
#                                   --gff3=on /home/skliver/data/genomes/fishing_cat/final.assembly.fasta > fishing_cat_augustus_predictions.gff
#
#cd ~/workdir/cats/leopard_cat
#~/soft/augustus-3.2.1/bin/augustus --species=Felis_catus --AUGUSTUS_CONFIG_PATH=/home/skliver/soft/augustus-3.2.1/config/ \
#                                   --strand=both --genemodel=complete \
#                                   --gff3=on /home/skliver/data/genomes/leopard_cat/final.assembly.fasta > leopard_cat_augustus_predictions.gff
#
#

cd ~/workdir/cats/fishing_cat/parallel_augustus
~/soft/MAVR/scripts/annotation/parallel_augustus.py -i /home/skliver/data/genomes/fishing_cat/final.assembly.fasta \
                                                    -o fishing_cat_augustus_cat_model.gff -t 10 -s Felis_catus \
                                                    -c /home/skliver/soft/augustus-3.2.1/config/

cd ~/workdir/cats/leopard_cat/parallel_augustus
~/soft/MAVR/scripts/annotation/parallel_augustus.py -i /home/skliver/data/genomes/leopard_cat/final.assembly.fasta \
                                                    -o leopard_cat_augustus_cat_model.gff -t 10 -s Felis_catus \
                                                    -c /home/skliver/soft/augustus-3.2.1/config/

cd ~/workdir/cats/caracal/parallel_augustus
~/soft/MAVR/scripts/annotation/parallel_augustus.py -i /home/skliver/data/genomes/caracal/final.assembly.fasta \
                                                    -o caracal_augustus_cat_model.gff -t 10 -s Felis_catus \
                                                    -c /home/skliver/soft/augustus-3.2.1/config/

