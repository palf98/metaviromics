#-*- coding: UTF-8 -*-
import sys
import re

#antes:
#for FILE in /home/paualico/metavirome/FINAL/*/*/*_1_P_Trimmomatic.fastq.gz; do ACC="${FILE/_1_P_Trimmomatic.fastq.gz/}"; dir=${FILE/_1_P_Trimmomatic.fastq.gz/_kallisto}; echo "$ACC" | awk -F "/" '{print $NF}'; grep TRINITY "$dir"/abundance.tsv | wc -l; done > ../metavirome/FINAL/ncontigs.txt



nline = 1
print('acc,ncontigs')
for line in open("/home/paualico/metavirome/FINAL/ncontigs.txt"):
    if (nline % 2) != 0:
        a = line.strip('\n')
    elif (nline % 2) == 0:
        b = line.strip('\n')
        print(a + ',' + b)
    nline += 1