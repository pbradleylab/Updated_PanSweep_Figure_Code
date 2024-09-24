#!/bin/bash
#Getting the nucleotide sequences for each gene

#Paths#
G_100060="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100060/centroids.ffn"
G_100078="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100078/centroids.ffn"
G_100271="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100271/centroids.ffn"
G_101380="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/101380/centroids.ffn"
G_100217="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100217/centroids.ffn"


Save_Path="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/"
SV_G_100060="New_Gene_sequences_100060.txt"
SV_G_100078="New_Gene_sequences_100078.txt"
SV_G_100271="New_Gene_sequences_100271.txt"
SV_G_101380="Gene_sequences_101380.txt"
SV_G_100217="Gene_sequences_100217.txt"

PT_G_100060="${Save_Path}${SV_G_100060}"
PT_G_100078="${Save_Path}${SV_G_100078}"
PT_G_100271="${Save_Path}${SV_G_100271}"
PT_G_101380="${Save_Path}${SV_G_101380}"
PT_G_100217="${Save_Path}${SV_G_100217}"

#######################################################

Genes_100060=("249173_01879" "148769_00831" "108348_02595" "210755_00357" "239171_01612" "108323_00300" "030659_00967" "001288_02675" "033855_01832" "152430_02393" "038364_02001" "235946_01311" "000216_02420" "057388_00211" "083468_01778" "230652_00809" "051001_01727" "033855_01831" "036212_01418" "052248_01872" "011927_00007")

for gene in "${Genes_100060[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100060" | sed '$d' >> "$PT_G_100060"
done
######################################################

Genes_100078=("109369_01803")


for gene in "${Genes_100078[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100078" | sed '$d' >> "$PT_G_100078"
done
#####################################################
Genes_100271=("000117_02863" "001770_01542" "155662_03708" "004804_02592" "198550_01569" "000355_03097" "027032_02679" "083479_02753" "137889_03520" "000117_00039" "006614_00504" "115855_01382" "155827_01231" "125279_02764" "000117_03451")


for gene in "${Genes_100271[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100271" | sed '$d' >> "$PT_G_100271"
done
#####################################################
Genes_101380=("000135_01093"   "005482_00240"   "187977_00187"   "076937_00800"   "169515_01624"  )


for gene in "${Genes_101380[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_101380" | sed '$d' >> "$PT_G_101380"
done
###################################################
Genes_100217=("047729_02814"   "031386_02238"   "046826_00950")


for gene in "${Genes_100217[@]}"; do
	sed -n "/^>UHGG${gene}/,/^>/p" "$G_100217" | sed '$d' >> "$PT_G_100217"
done
###################################################
