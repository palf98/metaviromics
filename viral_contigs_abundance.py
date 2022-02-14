# This script takes VF postprocessed csv files and VS2 tsv files as input and uses ID's to compute
# family abundance based on Kallisto's estimated counts and Diamond's taxonomic classification of contigs.
# Families are identified from taxIDs running taxonkit in lineage mode.

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
parser.add_argument("--vf", default=None, nargs = '+', help= "Absolute path to VirFinder post-processed output, i.e. a csv file containing a selection of contigs based on p-value or corrected p-value. One or more files can be provided one after the other (--vf X Y) or using bash wildcards.")
parser.add_argument("--vs", default=None, help= "Absolute path to VirSorter2 output (final-viral-combined.fa).")
parser.add_argument("--mode", default = "viral", help = "Mode of handling Diamond output: ignore contigs without viral hits ('--mode viral', default) or take full table ('--mode full', can be VERY SLOW)")
# action = store_true: if --verbose is called, args.verbose is True. Else it's False.
parser.add_argument("--verbose", action = "store_true", help= "Verbose functionality")
#output argument?

args = parser.parse_args()

if args.verbose:
	print('Run Started in '+ str(args.mode) +' mode: ' + time.ctime(time.time()))


#### PREVIOUS, COMMON STEP: select viral rows and add viral family to diamond tsv ####

if os.path.isfile(str(args.diamond)[0:-4] + '_families.tsv') == False:
    # get contig abundance from kallisto output file
    contig2klst = {}
    for klstline in open(args.kallisto):
    	klstline = klstline.strip('\n')
    	contig2klst[klstline.split('\t')[0]] = '\t' + klstline.split('\t')[3] + '\t' + klstline.split('\t')[4]

    # get full taxonomy lineage w/ taxonkit only for viral hits
    ps = subprocess.Popen(('grep', 'Viruses', args.diamond), stdout=subprocess.PIPE)
    output = subprocess.check_output((args.taxonkit, "lineage", "--taxid-field", "9", "--data-dir", args.dmp_dir), stdin=ps.stdout)
    ps.wait()

    if args.verbose:
        print(time.ctime(time.time()) + ': Lineages recovered for viral contigs in ' + args.diamond)
        # convert class bytes output to string
        taxdata = output.decode('utf-8')

    tsv_full = open(str(args.diamond)[0:-4] + '_families.tsv','w')

    for line in open(args.diamond):
        line = line.strip('\n')
        if line.split('\t')[10] == "Viruses":
            for taxline in taxdata.split('\n'):
                # connect lines with same contig ID in diamond and taxonkit output and also by hit accession, as one contig can be aligned with different sequences
                if (line.split('\t')[0] == taxline.split('\t')[0]) and (line.split('\t')[1] == taxline.split('\t')[1]):
                    # ignore viral rows with incomplete taxonomic data
                    for tax in taxline.split('\t')[13].split(';'):
                    	if 'viridae' in tax:
                    		line += '\t' + tax
                    		break
                    #if (len(taxline.split('\t')[13].split(';')) > 6) and ('viridae' in taxline.split('\t')[13]):
                    #   line += '\t' + taxline.split('\t')[13].split(';')[6]
                    # for viral contigs of unknown family
                    else:
                        line += '\tUnclassified viruses'
                # NA for non-viral contigs
        else:
            line += '\tNA'

        tsv_full.write(line + contig2klst[line.split('\t')[0]] + '\n')

    if args.verbose:
        print(time.ctime(time.time()) + ': Diamond tsv output with viral family column written to ' + str(args.diamond)[0:-4] + '_families.tsv')
else:
    if args.verbose:
        print(str(args.diamond)[0:-4] + '_families.tsv' + ' already exists. Using it for following steps.')


#### Generate tsv with contigs with any viral hit called by Diamond (viral mode) ####

if args.mode == "viral":
	viral_tsv = open(str(args.diamond)[0:-4] + '_families_contigs_viral_hit.tsv', 'w')
	# contig stored if any aligned sequence to the contig is viral
	contig_is_viral = set()
	for line in open(str(args.diamond)[0:-4] + '_families.tsv'):
		line = line.strip('\n')
		if line.split('\t')[0] not in contig_is_viral:
			if line.split('\t')[10] == "Viruses":
				contig_is_viral.add(line.split('\t')[0])
				
	for line in open(str(args.diamond)[0:-4] + '_families.tsv'):
		line = line.strip('\n')
		if (line.split('\t')[0] in contig_is_viral) and (line.split('\t')[10] == 'Viruses'):
			viral_tsv.write(line + '\n')


#### Given a selection of contigs (in a python set), sum est_counts and tpm per family ####
def compute_abundances(contigs, setname, mode):
	contig2family = {} # value: family with lowest e-value among contig's hits
	contig2est_counts = {}
	contig2tpm = {}
	if args.mode == 'full':
		if args.verbose:
			print(time.ctime(time.time()) + ': Looking for '+ str(args.vs) +' contigs in ' + str(args.diamond)[0:-4] + '_families.tsv' + '. This might take a while.')
		for contig in contigs:
			lowesteval = 1000
			for line in open(str(args.diamond)[0:-4] + '_families.tsv'):
				line = line.strip('\n')
				if line.split('\t')[0] == contig:
					if contig not in contig2tpm:
						lowesteval = float(line.split('\t')[6])
						contig2est_counts[contig] = line.split('\t')[14]
						contig2tpm[contig] = line.split('\t')[15]
						contig2family[contig] = line.split('\t')[13]
					#get family of lowest e-val score hit per contig
					else:
						if float(line.split('\t')[6]) < lowesteval: # contig hits are ordered by default from lowest e-value to highest, this is not necessary if input is a normal diamond file
							lowesteval = float(line.split('\t')[6])
							contig2family[contig] = line.split('\t')[13]

	elif args.mode == 'viral':
		if args.verbose:
			print(time.ctime(time.time()) + ': Looking for '+ str(setname) +' contigs in ' + str(args.diamond)[0:-4] + '_families_contigs_viral_hit.tsv.')
		for contig in contigs:
			lowesteval = 1000
			for line in open(str(args.diamond)[0:-4] + '_families_contigs_viral_hit.tsv'):
				line = line.strip('\n')
				if (line.split('\t')[0] == contig) and (line.split('\t')[10] == 'Viruses'):
					if contig not in contig2tpm:
						lowesteval = float(line.split('\t')[6])
						contig2est_counts[contig] = line.split('\t')[14]
						contig2tpm[contig] = line.split('\t')[15]
						contig2family[contig] = line.split('\t')[13]
					#get family of lowest e-val score hit per contig
					else:
						if float(line.split('\t')[6]) < lowesteval: # contig hits are ordered by default from lowest e-value to highest, this is not necessary if input is a normal diamond file
							lowesteval = float(line.split('\t')[6])
							contig2family[contig] = line.split('\t')[13]

		fam2est_counts = {}
		fam2tpm = {}
		for contig in contig2family:
			if contig2family[contig] not in fam2est_counts:
				fam2est_counts[contig2family[contig]] = float(contig2est_counts[contig])
				fam2tpm[contig2family[contig]] = float(contig2tpm[contig])
			else:
				fam2est_counts[contig2family[contig]] += float(contig2est_counts[contig])
				fam2tpm[contig2family[contig]] += float(contig2tpm[contig])

		fam_out = open(str(setname) + '_family_abundance.tsv','w')
		for family in sorted(fam2est_counts.keys()):
			fam_out.write(family + '\t' + str(fam2est_counts[family]) + '\t' + str(fam2tpm[family]) + '\n' )
		if args.verbose:
			print(time.ctime(time.time()) + ': Family abundances written to ' + setname + '_family_abundance.tsv')
		fam_out.close()



#### VirFinder contig set(s) ####

if args.vf:
	filenames2contigset = {}
	for vffile in args.vf:
		#name = '.'.join(os.path.split(vffile)[-1].split('.')[:-1]) #strip file format as last field splitting by '.'
		name = '.'.join(vffile.split('.')[:-1])
		filenames2contigset[name] = set()
		for line in open(vffile):
			line = line.strip('\n')
			if line[0] == '>':
				line = line.strip('>')
				filenames2contigset[name].add(line.split(' ')[0])
			else:
				continue
		if args.verbose:
			print(time.ctime(time.time()) + ': Contigs recovered for ' + vffile)
		
		compute_abundances(filenames2contigset[name], name, args.mode)	



#### VirSorter2 contig set ####

if args.vs:
	vscontigs = {}
	vscontigs[args.vs] = set()
	for line in open(args.vs):
		line = line.strip('\n')
		if line[0] == '>':
			line = line.strip('>')
			vscontigs[args.vs].add(line.split('||')[0])
		else:
			continue
	if args.verbose:
		print(time.ctime(time.time()) + ': Contigs recovered for ' + args.vs)

	compute_abundances(vscontigs[args.vs], '.'.join(args.vs.split('.')[:-1]), args.mode)



#### Abundances straight out from Diamond output ####

if args.mode == 'viral':
	compute_abundances(contig_is_viral, str(args.diamond)[0:-4], args.mode)
