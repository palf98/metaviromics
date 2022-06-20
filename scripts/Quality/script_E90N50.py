# This script needs: 1) fasta assembly and 2) transcript abundance matrix (kallisto's abundance.tsv)

import sys
import subprocess
import argparse
import os

parser = argparse.ArgumentParser(description='A script to extract ExN50 statistic and relate it to run accession ID')

parser.add_argument("--assembly", required = "True", default = None, help = "Path to fasta assembly")
parser.add_argument("--abundance", required = "True", default = None, help = "Path to transcript abundance table (e.g. kallisto's 'abundance.tsv')")
parser.add_argument("--x", default = 90, help = "N50 will be estimated for top x percentage most abundant transcripts, 90 by default. Integers between 0 and 100.")
parser.add_argument("--trinity_path", required = "True", default = None, help = "Path to Trinity bin directory")
args = parser.parse_args()


# get accession ( substring RR has to be absent in the path besides from run ID)
acc = args.assembly[(args.assembly.find("RR")) - 1:][:args.assembly[(args.assembly.find("RR")) - 1:].find("_")]

ps = str(subprocess.check_output((str(args.trinity_path+'contig_ExN50_statistic.pl'), args.abundance, args.assembly)), 'utf-8')

for line in ps.split('\n'):
    line = line.strip('\n')
    if line.split('\t')[0] in ['Ex', '']:
        continue
    elif int(line.split('\t')[0]) == int(args.x):
        print(acc + '\t' + line)







