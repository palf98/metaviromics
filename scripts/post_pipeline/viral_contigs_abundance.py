### EN PRIMER LUGAR

# Se han quitado las comparaciones que hacia el script original para comprobar si el hit era viral, 
# ahora asumimos que todos son virales al haber restringido 'nr' por el TaxID 'Viruses'.

# Se han modificado los indices de las columnas para encajar con 
#'qseqid sseqid qstart qend sstart send evalue bitscore length pident staxids sscinames sskingdoms skingdoms sphylums'

# El output se genera en la misma carpeta que los inputs, es decir, en PHYLUM/SPECIES/.

# Run with Python 3

import sys
import subprocess
import argparse
import os
import time

parser = argparse.ArgumentParser(description='A script to combine Diamond, VirFinder, VirSorter2 and Kallisto outputs to obtain viral family abundances in metatranscriptomics! For each contig, viral family of the viral hit with lowest e-value and Kallisto abundance metrics are stored. Afterwards, metrics are summed per family to obtain an overview of viral diversity in the sample.')

parser.add_argument("--diamond", required = "True", default = None, help = "Absolute path to Diamond tsv output generated with 'qseqid sseqid qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums' arguments. REQUIRED.")
parser.add_argument("--kallisto", required = "True", default = None, help = "Absolute path to Kallisto abundance.tsv file. REQUIRED.")
parser.add_argument("--taxonkit", required = "True", default = None, help = "Absolute path to Taxonkit binary. REQUIRED.")
parser.add_argument("--dmp_dir", required = "True", default = None, help = "Absolute path to directory containing .dmp files for the database used in Diamond. REQUIRED.")
# action = store_true: if --verbose is called, args.verbose is True. Else it's False.
parser.add_argument("--verbose", action = "store_true", help= "Verbose functionality")


args = parser.parse_args()

#### PREVIOUS, COMMON STEP: select viral rows and add viral family to diamond tsv ####

if os.path.isfile(str(args.diamond)[0:-4] + '_families.tsv') == False:
    # get contig abundance from kallisto output file
    contig2klst = {}
    for klstline in open(args.kallisto):
        klstline = klstline.strip('\n')
        contig2klst[klstline.split('\t')[0]] = '\t' + klstline.split('\t')[3] + '\t' + klstline.split('\t')[4]

    # get full taxonomy lineage w/ taxonkit only for viral hits
    ps = subprocess.Popen(('grep', 'Viruses', args.diamond), stdout=subprocess.PIPE)
    output = subprocess.check_output((args.taxonkit, "lineage", "--taxid-field", "11", "--data-dir", args.dmp_dir), stdin=ps.stdout)
    ps.wait()

    if args.verbose:
        print(time.ctime(time.time()) + ': Lineages recovered for viral contigs in ' + args.diamond)
        # convert class bytes output to string
        taxdata = output.decode('utf-8')

    tsv_full = open(str(args.diamond)[0:-4] + '_families.tsv','w')

    for taxline in taxdata.split('\n'):
        if len(taxline.split('\t')) < 10:
            break
        for tax in taxline.split('\t')[15].split(';'):
            if 'viridae' in tax:
                taxline += '\t' + tax
                break
        else:
            taxline += '\tUnclassified viruses'
        tsv_full.write(taxline + contig2klst[taxline.split('\t')[0]] + '\n')

    if args.verbose:
    	print(time.ctime(time.time()) + ': Diamond tsv output with viral family column written to ' + str(args.diamond)[0:-4] + '_families.tsv')
else:
    if args.verbose:
        print(str(args.diamond)[0:-4] + '_families.tsv' + ' already exists. Using it for following steps.')

#### Given a selection of contigs (in a python set), sum est_counts and tpm per family ####
def compute_abundances(setname):
	contig2family = {} # value: family with lowest e-value among contig's hits
	contig2est_counts = {}
	contig2tpm = {}
	if args.verbose:
		print(time.ctime(time.time()) + ': Looking for '+ str(setname) +' contigs in ' + str(args.diamond)[0:-4] + '_families.tsv' + '. This might take a while.')
	for line in open(str(args.diamond)[0:-4] + '_families.tsv'):
		contig = line.split('\t')[0]
		lowesteval = 1000
		line = line.strip('\n')
		if contig not in contig2tpm:
			lowesteval = float(line.split('\t')[6])
			contig2est_counts[contig] = line.split('\t')[17]
			contig2tpm[contig] = line.split('\t')[18]
			contig2family[contig] = line.split('\t')[16]
		#get family of lowest e-val score hit per contig
		else:
			if float(line.split('\t')[6]) < lowesteval: # contig hits are ordered by default from lowest e-value to highest, this is not necessary if input is a normal diamond file
				lowesteval = float(line.split('\t')[6])
				contig2family[contig] = line.split('\t')[16]

	fam2est_counts = {}
	fam2tpm = {}
	for contig in contig2family:
		if contig2family[contig] not in fam2est_counts:
			fam2est_counts[contig2family[contig]] = float(contig2est_counts[contig])
			fam2tpm[contig2family[contig]] = float(contig2tpm[contig])
		else:
			fam2est_counts[contig2family[contig]] += float(contig2est_counts[contig])
			fam2tpm[contig2family[contig]] += float(contig2tpm[contig])
	# SPLIT OUTPUT TYPE
	fam_out1 = open(str(setname) + '_family_abundance_est_counts.tsv','w')
	fam_out2 = open(str(setname) + '_family_abundance_tpm.tsv','w')
	for i in range(3):
		if i == 0:
			headerline = ''
			for family in sorted(fam2est_counts.keys()):
				headerline += family + '\t'
			fam_out1.write(headerline.strip('\t') + '\n')
			fam_out2.write(headerline.strip('\t') + '\n')
		if i == 1:
			est_counts_line = ''
			for family in sorted(fam2est_counts.keys()):
				est_counts_line += str(fam2est_counts[family]) + '\t'
			fam_out1.write(est_counts_line.strip('\t') + '\n')
		if i == 2:
			tpm_line = ''
			for family in sorted(fam2est_counts.keys()):
				tpm_line += str(fam2tpm[family]) + '\t'
			fam_out2.write(tpm_line.strip('\t') + '\n')
	fam_out1.close()
	fam_out2.close()
	if args.verbose:
		print(time.ctime(time.time()) + ': Family abundances (est_counts) written to ' + setname + '_family_abundance_est_counts.tsv')
		print(time.ctime(time.time()) + ': Family abundances (tpm) written to ' + setname + '_family_abundance_tpm.tsv')



#### Abundances straight out from Diamond output ####

compute_abundances(str(args.diamond)[0:-4])

