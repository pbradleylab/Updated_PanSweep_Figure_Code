#!/bin/bash
#Getting the nucleotide sequences for each gene

#Paths#
G_100060="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100060/centroids.ffn"
G_100078="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100078/centroids.ffn"
G_100271="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100271/centroids.ffn"
G_102528="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/102528/centroids.ffn"


Save_Path="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/"
SV_G_100060="Gene_sequences_100060.txt"
SV_G_100078="Gene_sequences_100078.txt"
SV_G_100271="Gene_sequences_100271.txt"
SV_G_102528="Gene_sequences_102528.txt"

PT_G_100060="${Save_Path}${SV_G_100060}"
PT_G_100078="${Save_Path}${SV_G_100078}"
PT_G_100271="${Save_Path}${SV_G_100271}"
PT_G_102528="${Save_Path}${SV_G_102528}"

#######################################################

Genes_100060=("047117_02378"   "047117_02380"   "047117_02376"   "047117_02377"   "041746_02389"   "000216_02074"   "055467_01854"   "047117_02382"   "239171_01612"   "001288_02675"   "063307_00097"   "047117_02379"   "210793_02030"   "041746_02390"   "029873_02068"   "210928_01151")

for gene in "${Genes_100060[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100060" | sed '$d' >> "$PT_G_100060"
done
######################################################

Genes_100078=("064037_00369"   "154988_02076"   "016431_01839"   "108644_01598"   "188302_00118"   "151972_01002"   "195023_01285")


for gene in "${Genes_100078[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100078" | sed '$d' >> "$PT_G_100078"
done
#####################################################
Genes_100271=("200056_01835"   "032492_02428"   "200056_01834"   "200056_01833"   "200056_01832"   "155662_03633"   "155662_03778")


for gene in "${Genes_100271[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100271" | sed '$d' >> "$PT_G_100271"
done
###################################################
Gene_102528=("064838_00419")

for gene in "${Gene_102528[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_102528" | sed '$d' >> "$PT_G_102528"
done
