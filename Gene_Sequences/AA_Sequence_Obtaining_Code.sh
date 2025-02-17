#!/bin/bash


#SBATCH --job-name=AA_UHGP_100_Sequences
#SBATCH --time=12:00:00
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=1 
#SBATCH --mail-type=ALL
 


#Code to get aa sequence for each gene 
#Paths:
To_UHGP_100="/fs/project/bradley.720/db/uhgp/200624/uhgp-100"
Save_location="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/aa_sequences_PanSweep_Analysis.txt"

cd "$To_UHGP_100"

grep -A1 -e "041746_02389"   -e "156581_00568"   -e "246241_00488"   -e "200056_01834"   -e "186001_01421"   -e "016431_01839"   -e "149870_00780"   -e "210811_00453"   -e "181581_00481"   -e "210531_01353"   -e "258864_02598"   -e "200056_01833"   -e "041746_02388"   -e "151055_01601"   -e "055467_01854"   -e "004689_00515"   -e "058324_00427"   -e "030231_00195"   -e "152024_01664"   -e "063154_01644"   -e "060538_00391"   -e "017453_01069"   -e "000216_02074"   -e "155662_03633"   -e "099736_00986"   -e "007643_02005"   -e "017126_00406"   -e "192308_01194"   -e "010616_00971"   -e "154577_00206"   -e "063941_02093"   -e "029873_02068"   -e "210928_01151"   -e "063495_00342"   -e "200056_01835"   -e "064838_00419"   -e "047117_02377"   -e "047117_02380"   -e "032185_00468"   -e "108644_01598"   -e "176730_00858"  uhgp-100.faa >> "$Save_location"

