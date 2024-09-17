#!/bin/bash

#Code to get aa sequence for each gene 
#Paths:
To_UHGP_100="/fs/project/bradley.720/db/uhgp/200624/uhgp-100"
Save_location="/home/majernik.14/Updated_PanSweep_Paper_Figure_Code/aa_sequences_PanSweep_Analysis.txt"

cd "$To_UHGP_100"

grep -A1  uhgp-100.faa >> "$Save_location"

