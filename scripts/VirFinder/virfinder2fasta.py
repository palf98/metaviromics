# modified 17/02: now it works no matter if there is a q-value column as fdr005 and percentile30 logic vector columns are searched by column name and not index


import os
import sys

# absolute path to a single virfinder postprocessed csv
incsv = str(sys.argv[1])

infasta = '/' + '/'.join(incsv.split('/')[1:-2]) + '/' + incsv.split('/')[len(incsv.split('/')) - 1].split('_')[0] + '_trinity_out.Trinity.fasta'

selec30 = []
selecfdr = []

seqs = {}

output_fdr = open(str(incsv.split('_postprocess.csv')[0] + '_fdrless005.fasta'), 'w')
output_30 = open(str(incsv.split('_postprocess.csv')[0] + '_percentile30.fasta'), 'w')


# getting selected sequences following fdr or percentile30 criteria
firstline = 1
for line in open(incsv):
    line = line.strip('\n')
# ADDED
#
    if firstline == 1:
        for i in range(len(line.split(','))):
            if line.split(',')[i] == '"fdr.less0.05"':
                fdrfield = i
            if line.split(',')[i] == '"percentil30"':
                perc30field = i
        firstline = 0
#
    else:
        if line.split(',')[fdrfield] == 'TRUE':
            selecfdr.append(line.split(',')[0].strip('"'))    
            seqs[line.split(',')[0].strip('"')] = ''
        if line.split(',')[perc30field] == 'TRUE':
            selec30.append(line.split(',')[0].strip('"'))
            seqs[line.split(',')[0].strip('"')] = ''
        else:
            continue
        
firsttime = 1
# searching for selected sequences in Trinity.fasta file
for line in open(infasta):
    line = line.strip('\n')
    if line[0] == '>':
        if firsttime == 0:
            if fdr_coincidence:
                output_fdr.write('>' + seqname + '\n' + seqs[seqname] + '\n')
            if selec30_coincidence:
                output_30.write('>' + seqname + '\n' + seqs[seqname] + '\n')
        if firsttime:
            firsttime = 0
        seqname = line[1:len(line)]
        if seqname in selecfdr and seqname in selec30:
            fdr_coincidence = 1
            selec30_coincidence = 1
        elif seqname in selecfdr:
            fdr_coincidence = 1
            selec30_coincidence = 0
        elif seqname in selec30:
            fdr_coincidence = 0
            selec30_coincidence = 1
        else: 
            fdr_coincidence = 0
            selec30_coincidence = 0
    else:
        if (fdr_coincidence == 1) or (selec30_coincidence == 1):
            seqs[seqname] += line

#last sequence
if fdr_coincidence:
    output_fdr.write('>' + seqname + '\n' + seqs[seqname] + '\n')
if selec30_coincidence:
    output_30.write('>' + seqname + '\n' + seqs[seqname] + '\n')

output_fdr.close()
output_30.close()


#print(seqs)        
