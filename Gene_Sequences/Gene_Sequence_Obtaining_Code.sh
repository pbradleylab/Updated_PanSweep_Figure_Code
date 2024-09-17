#!/bin/bash
#Getting the nucleotide sequences for each gene

#Paths#
G_100060="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100060"
G_100078="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100078"
G_100271="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/100271"
G_102528="/fs/project/bradley.720/db/midasdb_uhgg/pangenomes/102528"

Save_Path = #Save path
Save_Path="/home/majernik.14/Update_PanSweep_Analysis/Updated_PanSweep_Figure_Code/Gene_Sequences/"
SV_G_100060="Gene_sequences_100060.txt"
SV_G_100078="Gene_sequences_100078.txt"
SV_G_100271="Gene_sequences_100271.txt"
SV_G_102528="Gene_sequences_102528.txt"

PT_G_100060="${Save_Path}${SV_G_100060}"
PT_G_100078="${Save_Path}${SV_G_100078}"
PT_G_100271="${Save_Path}${SV_G_100271}"
PT_G_102528="${Save_Path}${SV_G_102528}"


cd "$G_100060"

sed -n '/047117_02378/,/p/; /047117_02380/,/p/; /047117_02376/,/p/; /047117_02377/,/p/; /041746_02389/,/p/; /000216_02074/,/p/; /055467_01854/,/p/; /047117_02382/,/p/; /239171_01612/,/p/; /001288_02675/,/p/; /063307_00097/,/p/; /047117_02379/,/p/; /210793_02030/,/p/; /041746_02390/,/p/; /029873_02068/,/p/; /210928_01151/,/p/' centroids.ffn >> "$PT_G_100060"

cd "$G_100078"

sed -n '/064037_00369/,/p/; /154988_02076/,/p/; /016431_01839/,/p/; /108644_01598/,/p/; /188302_00118/,/p/; /151972_01002/,/p/; /195023_01285/,/p/' centroids.ffn >> "$PT_G_100078"

cd "$G_100271"

sed -n '/200056_01835/,/p/; /032492_02428/,/p/; /200056_01834/,/p/; /200056_01833/,/p/; /200056_01832/,/p/; /155662_03633/,/p/; /155662_03778/,/p/' centroids.ffn >> "$PT_G_100271"

cd "$G_102528"

sed -n '/064838_00419/,/p/' centroids.ffn >> "$PT_G_102528"

