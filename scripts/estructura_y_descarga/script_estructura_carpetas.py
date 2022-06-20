# Creates (phylum <- species) folder structure

import sys
import os


# excel database (tsv)
infile = sys.argv[1]

# folder that will contain phylum folders, script mkdirs it if non-existant
output_folder = str(sys.argv[2])
if output_folder[-1] != '/':
    output_folder += '/'
if os.path.isdir(output_folder) != True:
    os.mkdir(output_folder)


# comma-separated run accessions for each species
especies_acc = {}
especies_phyl = {}
especies_unicas = []

header = 1

for line in open(infile):
    if header == 1:
        header = 0
        continue
    else:
        line = line.strip('\n')
        especiedelafila = line.split('\t')[12]
        if especiedelafila in especies_unicas:
            especies_acc[especiedelafila] += ',' + line.split('\t')[1]
        else:
            # si la especie no esta en unicas, la anyadimos
            especies_unicas.append(especiedelafila)
            # abrimos entrada en los diccionarios para las especies unicas
            especies_acc[especiedelafila] = line.split('\t')[1]
            especies_phyl[especiedelafila] = line.split('\t')[8]
            
especies_unicas = sorted(especies_unicas)               
phyl_set = set(sorted(especies_phyl.values()))
phyl_list = list(phyl_set)

# PHYLUM superfolder

for phylum in phyl_list:
    phylum = phylum.strip('"')
    os.mkdir(output_folder + phylum.upper())
    for especie in especies_unicas:
        phyl2 = especies_phyl[especie].strip('"')
        if phyl2 == phylum:
            especie = especie.strip('"')
            # nombres de archivos sin espacios, las especies no tienen '_' en su nombre, no habra inconsistencias
            nombre_especie = especie.replace(' ','_')
            os.mkdir(output_folder + phylum.upper() + '/' + nombre_especie)
            out = open(output_folder + phylum.upper() + '/' + nombre_especie + '/' + nombre_especie + '.acc', 'w')
            for acc in especies_acc[especie].split(','):
                acc = acc.strip('"')
                out.write(acc + '\n')
            out.close()
            