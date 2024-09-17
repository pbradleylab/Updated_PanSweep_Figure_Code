#!/bin/bash

#Paths:
To_UHGP_100="/fs/project/bradley.720/db/uhgp/200624/uhgp-100"
Save_location="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/Uhgp_100_to_Gene_Numbers.tsv"


#To obtain the UHGP-100 cluster_ids per gene 
cd "$To_UHGP_100"

grep -e "064037_00369" -e "047117_02378" -e "154988_02076" -e "016431_01839" -e "200056_01835" -e "108644_01598" -e "047117_02380" -e "032492_02428" -e "047117_02376" -e "047117_02377" -e "200056_01834" -e "041746_02389" -e "200056_01833" -e "000216_02074" -e "055467_01854" -e "047117_02382" -e "239171_01612" -e "188302_00118" -e "151972_01002" -e "001288_02675" -e "063307_00097" -e "200056_01832" -e "047117_02379" -e "195023_01285" -e "155662_03633" -e "155662_03778" -e "210793_02030" -e "041746_02390" -e "029873_02068" -e "210928_01151" -e "064838_00419" uhgp-100.tsv >> "$Save_location"
