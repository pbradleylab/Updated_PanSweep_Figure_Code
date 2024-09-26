#!/bin/bash

#Get presence/absence for each gene in each genomes that makes up a species

#Paths#
G_100060="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00060/pan-genome/genes_presence-absence_locus.csv"
G_100078="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00078/pan-genome/genes_presence-absence_locus.csv"
G_100271="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-002/MGYG-HGUT-00271/pan-genome/genes_presence-absence_locus.csv"
G_102528="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-025/MGYG-HGUT-02528/pan-genome/genes_presence-absence_locus.csv"
G_100217="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-002/MGYG-HGUT-00217/pan-genome/genes_presence-absence_locus.csv"
G_101380="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-013/MGYG-HGUT-01380/pan-genome/genes_presence-absence_locus.csv"

#Save path
Save_Path="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
SV_G_100060="Sig_Genes_Pangenomes_100060.csv"
SV_G_100078="Sig_Genes_Pangenomes_100078.csv"
SV_G_100271="Sig_Genes_Pangenomes_100271.csv"
SV_G_102528="Sig_Genes_Pangenomes_102528.csv"
SV_G_100217="Sig_Genes_Pangenomes_100217.csv"
SV_G_101380="Sig_Genes_Pangenomes_101380.csv"

PT_G_100060="${Save_Path}${SV_G_100060}"
PT_G_100078="${Save_Path}${SV_G_100078}"
PT_G_100271="${Save_Path}${SV_G_100271}"
PT_G_102528="${Save_Path}${SV_G_102528}"
PT_G_100217="${Save_Path}${SV_G_100217}"
PT_G_101380="${Save_Path}${SV_G_101380}"

#######################################################

Genes_100060=("249173_01879"   "047117_02378"   "228006_01743"   "148769_00831"   "047117_02380"   "148794_01632"   "047117_02376"   "047117_02377"   "041746_02389"   "108348_02595"   "192308_01194"   "000216_02074"   "210755_00357"   "152466_01649"   "055467_01854"   "047117_02383"   "047117_02382"   "181581_00481"   "239171_01612"   "032185_00468"   "108323_00300"   "030659_00967"   "001288_02675"   "041746_02388"   "033855_01832"   "063307_00097"   "047117_02379"   "152430_02393"   "038364_02001"   "210793_02030"   "235946_01311"   "041746_02390"   "047117_02375"   "047117_02381"   "000216_02420"   "057388_00211"   "083468_01778"   "230652_00809"   "051001_01727"   "033855_01831"   "036212_01418"   "052248_01872"   "258864_02598"   "029873_02068"   "210928_01151"   "011927_00007" )

H060=$(head -n 1 $G_100060)  
echo "\"UHGG_Gene\",$H060" > "$PT_G_100060"

for gene in "${Genes_100060[@]}"; do
	G_060=$(grep "$gene" $G_100060)
        echo "\"$gene\",$G_060">> "$PT_G_100060"
done
######################################################

Genes_100078=("064037_00369"   "154988_02076"   "016431_01839"   "108644_01598"   "109369_01803"   "063379_01966"   "188302_00118"   "151972_01002"   "195023_01285")

H271=$(head -n 1 $G_100078)

echo "\"UHGG_Gene\",$H271" > "$PT_G_100078"

for gene in "${Genes_100078[@]}"; do
	G_271=$(grep "$gene" $G_100078)
        echo "\"$gene\",$G_271" >> "$PT_G_100078"
done
#####################################################
Genes_100271=("000117_02863"   "001770_01542"   "200056_01835"   "032492_02428"   "155662_03708"   "200056_01834"   "200056_01833"   "004804_02592"   "198550_01569"   "000355_03097"   "027032_02679"   "083479_02753"   "137889_03520"   "000117_00039"   "200056_01832"   "155662_03633"   "155662_03778"   "006614_00504"   "115855_01382"   "155827_01231"   "125279_02764"   "000117_03451")

H694=$(head -n 1 $G_100271) 
echo "\"UHGG_Gene\",$H694" > "$PT_G_100271"

for gene in "${Genes_100271[@]}"; do
	G_694=$(grep "$gene" $G_100271) 
        echo "\"$gene\",$G_694" >> "$PT_G_100271"
done
###################################################
Gene_102528=("064838_00419")

H580=$(head -n 1 $G_102528)
echo "\"UHGG_Gene\",$H580" > "$PT_G_102528"

for gene in "${Gene_102528[@]}"; do
	G_580=$(grep "$gene" $G_102528) 
        echo "\"$gene\",$G_580" >> "$PT_G_102528"
done
#####################################################
Genes_101380=("000135_01093"   "005482_00240"   "187977_00187"   "076937_00800"   "169515_01624" )

H694=$(head -n 1 $G_101380) 
echo "\"UHGG_Gene\",$H694" > "$PT_G_101380"

for gene in "${Genes_101380[@]}"; do
	G_694=$(grep "$gene" $G_101380) 
        echo "\"$gene\",$G_694" >> "$PT_G_101380"
done
###################################################
Genes_100217=("047729_02814"   "031386_02238"   "046826_00950" )

H694=$(head -n 1 $G_100217) 
echo "\"UHGG_Gene\",$H694" > "$PT_G_100217"

for gene in "${Genes_100217[@]}"; do
	G_694=$(grep "$gene" $G_100217) 
        echo "\"$gene\",$G_694" >> "$PT_G_100217"
done
###################################################

