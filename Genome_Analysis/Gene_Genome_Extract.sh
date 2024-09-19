#!/bin/bash

#Get presence/absence for each gene in each genomes that makes up a species

#Paths#
G_100060="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00060/pan-genome/genes_presence-absence_locus.csv"
G_100078="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00078/pan-genome/genes_presence-absence_locus.csv"
G_100271="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-002/MGYG-HGUT-00271/pan-genome/genes_presence-absence_locus.csv"
G_102528="/fs/project/bradley.720/db/uhgg/uhgg_catalogue/MGYG-HGUT-025/MGYG-HGUT-02528/pan-genome/genes_presence-absence_locus.csv"

#Save path
Save_Path="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Genome_Analysis/"
SV_G_100060="Sig_Genes_Pangenomes_100060.csv"
SV_G_100078="Sig_Genes_Pangenomes_100078.csv"
SV_G_100271="Sig_Genes_Pangenomes_100271.csv"
SV_G_102528="Sig_Genes_Pangenomes_102528.csv"

PT_G_100060="${Save_Path}${SV_G_100060}"
PT_G_100078="${Save_Path}${SV_G_100078}"
PT_G_100271="${Save_Path}${SV_G_100271}"
PT_G_102528="${Save_Path}${SV_G_102528}"

#######################################################

Genes_100060=("047117_02378"   "228006_01743"   "047117_02380"   "148794_01632"   "047117_02376"   "041746_02389"   "047117_02377"   "192308_01194"   "000216_02074"   "152466_01649"   "055467_01854"   "047117_02383"   "047117_02382"   "181581_00481"   "032185_00468"   "041746_02388"   "063307_00097"   "047117_02379"   "210793_02030"   "041746_02390"   "047117_02375"   "047117_02381"   "258864_02598"   "029873_02068"   "210928_01151" )

H060=$(head -n 1 $G_100060)  
echo "\"UHGG_Gene\",$H060" > "$PT_G_100060"

for gene in "${Genes_100060[@]}"; do
	G_060=$(grep "$gene" $G_100060)
        echo "\"$gene\",$G_060">> "$PT_G_100060"
done
######################################################

Genes_100078=("108644_01598"   "064037_00369"   "154988_02076"   "016431_01839"   "063379_01966"   "188302_00118"   "151972_01002"   "195023_01285"  )

H271=$(head -n 1 $G_100078)

echo "\"UHGG_Gene\",$H271" > "$PT_G_100078"

for gene in "${Genes_100078[@]}"; do
	G_271=$(grep "$gene" $G_100078)
        echo "\"$gene\",$G_271" >> "$PT_G_100078"
done
#####################################################
Genes_100271=("200056_01835"   "032492_02428"   "200056_01833"   "200056_01834"   "200056_01832"   "155662_03633"   "155662_03778")

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


