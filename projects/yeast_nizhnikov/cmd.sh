#!/usr/bin/env bash

cd /home/mahajrod/Dropbox/yeast/30.05.2017

~/Genetics/MAVR/scripts/file/remove_by_column_value.py -i all.samples.all.tab -c 0 -f mobile_elements_repeats_pseudogenes.ids  -o all.samples.all.no_mobile.tab

replace_column_value_by_syn.py -i all.samples.all.no_mobile.tab -c 0 -s /home/winstorage/nizhnikov/genes.syn -o all.samples.all.no_mobile.common_names.tab