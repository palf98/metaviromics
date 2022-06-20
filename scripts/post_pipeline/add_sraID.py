### EN SEGUNDO LUGAR

# Simplemente dar como argumento la carpeta que esta por encima de los filos,
# ahora solo guardaremos el identificador de SRA, todas las demás variables 
# las sacaremos del excel conectando con este identificador.


import sys
import os


path_over_phyla = sys.argv[1]

for phylum in os.listdir(path_over_phyla):
	if os.path.isfile(os.path.join(path_over_phyla, phylum)) == True:
		continue
	for sp in os.listdir(os.path.join(path_over_phyla, phylum)):
		if os.path.isfile(os.path.join(path_over_phyla, phylum, sp)) == True:
			continue
		for item in os.listdir(os.path.join(path_over_phyla, phylum, sp)):

	#### Diamond ####

			if 'diamond_f_family_abundance_tpm.tsv' in item:
				acc = item.split('_')[0]
				nline = 0
				line1 = ''
				line2 = ''
				for line in open(os.path.join(path_over_phyla, phylum, sp, item)):
					if nline == 0:
						line1 = line
						nline = 1
					else:
						line2 = line
				if len(line.split('\t')) < 2:
                                       output = open(os.path.join(path_over_phyla, phylum, sp, str(item.split('.')[0]) + '_accession_ID.tsv'), 'w')
                                       output.write('accession' + '\n' + acc + '\n')
                                       output.close()
                                else:
                                       output = open(os.path.join(path_over_phyla, phylum, sp, str(item.split('.')[0]) + '_accession_ID.tsv'), 'w')
                                       output.write('accession\t' + line1 +  acc + '\t' + line2)
                                       output.close()

				#if len(line.split('\t')) < 2:
				#	output = open(os.path.join(path_over_phyla, phylum, sp, str(item.split('.')[0]) + '_host_tax.tsv'), 'w')
				#	output.write('phylum\tspecies\taccession' + '\n' + phylum + '\t' + sp + '\t' + acc + '\n')
				#	output.close()
				#else:
				#	output = open(os.path.join(path_over_phyla, phylum, sp, str(item.split('.')[0]) + '_host_tax.tsv'), 'w')
				#	output.write('phylum\tspecies\taccession\t' + line1 + phylum + '\t' + sp + '\t' + acc + '\t' + line2)
				#	output.close()

