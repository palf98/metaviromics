# -*- coding: utf-8 -*-
### script preparado para un output de kaiju-addTaxonNames con todos los posibles argumentos -r 
### (superkingdom, phylum, class, order, family, genus, species)

# El objetivo es generar un solo .csv para todos los kaiju_tax.txt contenidos bajo un directorio, relacionando el perfil
# de abundancias con la accesion, la especie y el phylum. Ademas, se anyade la columna 'nreads' con el numero de lecturas 
# de cada metatranscriptoma, solo para el fichero con contajes (absabun).


import sys
import os
import subprocess

abs_path_to_superfolder = str(sys.argv[1])
if abs_path_to_superfolder[-1] != "/":
	abs_path_to_superfolder += "/"

final_fichero = sys.argv[2] # lo que va despues de la accession en el nombre, por ejemplo "_1_P_kaiju.txt"

# species for each phylum
sp2phyl = {}
# SRR accessions for each species name
srr2sp = {}
fams_global = []

# relaciona cada accesion con el diccionario de frecuencias de familias generado en el loop
acc2fams = {}
acc2filas_virales = {}
acc2filas_totales = {}
# esta lista de accesiones tiene como unica utilidad escribir las filas del csv en el orden alfabetico de phylum y especie que tenemos en el ordenador
accs = []

phyla = sorted(os.listdir(abs_path_to_superfolder))
if "CAT" in phyla:
	phyla.remove("CAT")
for phylum in phyla:
	if os.path.isfile(abs_path_to_superfolder + phylum):
		continue	
	for sp in sorted(os.listdir(abs_path_to_superfolder + phylum)):
		# relacion especie -> phylum
		sp2phyl[sp] = phylum
		if os.path.isdir(abs_path_to_superfolder + phylum + "/" + sp) and len(sp) > 1:
			for file in os.listdir(abs_path_to_superfolder + phylum + "/" + sp + "/"):
				if (final_fichero) in file:	
					acc = str(file)[0:str(file).find('_')]
					accs.append(acc)
					# relacion accesion -> especie
					srr2sp[acc] = sp
					filas_virales = 0
					fams = {}
					# abrir el fichero kaiju_tax
					for line in open(abs_path_to_superfolder + phylum + "/" + sp + "/" + file):
						line = line.strip('\n')
						fields = line.split('\t')
						taxons = fields[7].split('; ')
						#eliminar elemento vacio tras '; ' (formato kaiju-addTaxonNames)
						taxons = taxons[:-1]
						# recuento de filas donde el superkingdom es 'Viruses', descartamos coincidencias eucariotas
						if taxons[0] == 'Viruses':
							filas_virales += 1
							# sumamos 1 a la frecuencia de la familia coincidente
							if taxons[4] not in fams.keys():
								fams[taxons[4]] = 1
								if taxons[4] not in fams_global:
									fams_global.append(taxons[4])
							else:
								fams[taxons[4]] += 1
					# frecuencia absoluta de cada familia de virus para la accesion
					acc2fams[acc] = fams
					# frecuencia absoluta total de secuencias clasificadas virales para la accesion
					acc2filas_virales[acc] = filas_virales
					# numero de lecturas (numero de lineas en acc_kaiju.txt)
					acc2filas_totales[acc] = subprocess.check_output('wc -l ' + str(abs_path_to_superfolder + phylum + "/" + sp + "/" + acc + final_fichero), shell = True).split(' ')[0]

### Escritura 
output = open(str(abs_path_to_superfolder) + "/relabun_kaiju.csv", "w")
output2 = open(str(abs_path_to_superfolder) + "/absabun_kaiju.csv", "w")
## Escritura de la linea HEADER (Phylum, Species, Accession, Virus family 1, Virus family 2, ...)
flag = 1
for entry in sorted(fams_global):
        if flag:
                #output.write("phylum,species,accession," + entry)
                output2.write("phylum,species,accession,nreads," + entry)
                flag = 0
        else:
                #output.write("," + entry)
                output2.write("," + entry)

#output.write("\n")
output2.write("\n")

# para cada accesion se escribe una fila, respetando el orden de columnas que dicta sorted(fams_global)
for acc in accs:
	flag = 1
	for family in sorted(fams_global):
		if family in acc2fams[acc].keys():
#			proporcion = float(acc2fams[acc][family]) / float(acc2filas_virales[acc])
			counts = acc2fams[acc][family]
		else:
#			proporcion = 0
			counts = 0
		if flag:
#			output.write(sp2phyl[srr2sp[acc]] + "," + srr2sp[acc] + "," + acc + "," + str(proporcion))
#			output2.write(sp2phyl[srr2sp[acc]] + "," + srr2sp[acc] + "," + acc + "," + str(counts))
			output2.write(sp2phyl[srr2sp[acc]] + "," + srr2sp[acc] + "," + acc + "," + str(acc2filas_totales[acc]) + "," + str(counts))
			flag = 0
		else:
#        		output.write("," + str(proporcion))
        		output2.write("," + str(counts))
#	output.write("\n")
	output2.write("\n")
	
#output.close()
output2.close()		