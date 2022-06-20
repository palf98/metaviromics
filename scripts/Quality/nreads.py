#-*- coding: UTF-8 -*-
import sys
import re


log = sys.argv[1]

ptrn1 = "clumpify start"
ptrn2 = "clumpify end"
ptrn3 = "Trimmomatic start"
ptrn4 = "Trimmomatic end"

# guardamos en acc_clump la parte del log de clump referente a una accession
acc_clump = {}
acc_trimmo = {}
acumulacion = ''
for line in open(log):
    if re.search(ptrn1, line):
        acc = line.split(' ')[0].split('/')[-1]
        acumulacion = line
    elif re.search(ptrn2, line):
        acc_clump[acc] = acumulacion + line
        acumulacion = ''
    else:
        acumulacion += line
    if re.search(ptrn3, line):
        acc = line.split(' ')[0].split('/')[-1]
        acumulacion = line
    elif re.search(ptrn4, line):
        acc_trimmo[acc] = acumulacion + line
    else:
        acumulacion += line
        
acc_clump_input_reads = {}
acc_clump_output_reads = {}
acc_trimmo_input_reads = {}
acc_trimmo_output_reads = {}

print('acc,clumpify input reads,clumpify output reads,Trimmomatic input reads,Trimmomatic output reads')

for acc in acc_clump.keys():
    for line in acc_clump[acc].split('\n'):
        if re.search("Reads In:", line):
            acc_clump_input_reads[acc] = line.split(' ')[-1]
        if re.search("Reads Out:", line): 
            acc_clump_output_reads[acc] = line.split(' ')[-1]
for acc in acc_trimmo.keys():
    for line in acc_trimmo[acc].split('\n'):
        if re.search("Both Surviving:", line):
            acc_trimmo_input_reads[acc] = line.split(': ')[1].split(' ')[0]
            acc_trimmo_output_reads[acc] = line.split(': ')[2].split(' ')[0]


for acc in acc_clump.keys():
    if acc in acc_clump_input_reads.keys():
        print(acc +','+acc_clump_input_reads[acc] + ',' + acc_clump_output_reads[acc] + ',' + 
        acc_trimmo_input_reads[acc] + ',' + acc_trimmo_output_reads[acc])
    else:
        print(acc + ',killed,killed,killed,killed')
        
        
        
        
        
        
        