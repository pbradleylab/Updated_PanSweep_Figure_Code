#!/bin/bash
#Getting the nucleotide sequences for each gene

cd /fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100060 

sed -n '/041746_02389/,/>/p; /047117_02376/,/>/p; /047117_02378/,/>/p; /047117_02379/,/>/p; /055467_01854/,/>/p; /152466_01649/,/>/p;  /210793_02030/,/>/p;  /047117_02377/,/>/p;  /029873_02068/,/>/p;  /210928_01151/,/>/p' centroids.ffn >> /home/majernik.14/PanSweep_Paper_Figure_Code/Gene_Sequence_code/PanSweep_Figure_Code/100060_gene_sequences.txt

cd /fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100271

sed -n '/200056_01832/,/>/p; /200056_01833/,/>/p; /200056_01834/,/>/p; /200056_01835/,/>/p' centroids.ffn >> /home/majernik.14/PanSweep_Paper_Figure_Code/Gene_Sequence_code/PanSweep_Figure_Code/100271_gene_sequences.txt

cd /fs/project/bradley.720/db/midasdb_uhgg/pangenomes/102580

sed -n '/008792_01807/,/>/p' centroids.ffn >> /home/majernik.14/PanSweep_Paper_Figure_Code/Gene_Sequence_code/PanSweep_Figure_Code/102580_gene_sequences.txt

cd /fs/project/bradley.720/db/midasdb_uhgg/pangenomes/103694

sed -n '/180526_01821/,/>/p' centroids.ffn >> /home/majernik.14/PanSweep_Paper_Figure_Code/Gene_Sequence_code/PanSweep_Figure_Code/103694_gene_sequences.txt

