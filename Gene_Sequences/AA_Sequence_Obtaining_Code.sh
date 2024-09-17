#!/bin/bash

#Code to get aa sequence for each gene 
#Paths:
To_UHGP_100="/fs/project/bradley.720/db/uhgp/200624/uhgp-100"
Save_location="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/aa_sequences_PanSweep_Analysis.txt"

cd "$To_UHGP_100"

grep -A1  uhgp-100.faa >> "$Save_location"

