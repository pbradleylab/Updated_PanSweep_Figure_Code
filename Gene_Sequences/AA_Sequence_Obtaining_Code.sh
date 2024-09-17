#!/bin/bash

#Code to get aa sequence for each gene 

cd /fs/project/bradley.720/db/uhgp/200624/uhgp-100

grep -A1  uhgp-100.faa >> /home/majernik.14/PanSweep_Paper_Figure_Code/Gene_Sequence_code/PanSweep_Figure_Code/aa_sequences_PanSweep_Analysis.txt

