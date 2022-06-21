#!/bin/bash
#
#SBATCH --job-name=------
#SBATCH --output=------
#
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=1000

####
#### CARPETAS ####
####

python /home/paualico/metavirome/scripts/script_structure.py /home/paualico/metavirome/FINAL_resto_filos.tsv /home/paualico/metavirome/FINAL/

####
#### DESCARGA ####
####

python /home/paualico/metavirome/scripts/script_ferqd.py /home/paualico/metavirome/FINAL/ /home/paualico/sratoolkit.2.11.2-ubuntu64/bin/ 25

####
#### PREPROCESADO DE READS (DEDUPLICACION, FILTRAR rRNA, REPARAR EMPAREJADO, CALIDAD) ####
#### 25 cores

module load anaconda/anaconda3

# Clumpify
source activate bbmap_conda

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1.fastq.gz; do FILE_2="${FILE/_1.fastq.gz/_2.fastq.gz}"; echo "${FILE/_1.fastq.gz/} clumpify start"; date; /home/paualico/.conda/envs/bbmap_conda/bin/clumpify.sh in="$FILE" in2="$FILE_2" out="${FILE/_1.fastq.gz/_fwd_dedup.fq.gz}" out2="${FILE_2/_2.fastq.gz/_rev_dedup.fq.gz}" dedupe=t groups=auto lowcomplexity=f -Xmx50g; echo "${FILE/_1.fastq.gz/} clumpify end"; date; done

conda deactivate


# SortMeRNA with --paired_in to avoid reads paired to COI or rRNA hits in '--other' output
source activate sortmerna

for FILE in /home/paualico/metavirome/FINAL/*/*/*_fwd_dedup.fq.gz; do FILE_2="${FILE/_fwd_dedup.fq.gz/_rev_dedup.fq.gz}"; echo "${FILE/_1.fastq.gz/} sortmerna start"; date; /home/paualico/.conda/envs/sortmerna/bin/sortmerna -ref /home/paualico/.conda/envs/sortmerna/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta -ref /home/paualico/.conda/envs/sortmerna/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta -ref /home/paualico/.conda/envs/sortmerna/sortmerna/data/rRNA_databases/COI.fasta -reads "$FILE" -reads "$FILE_2" -threads 25 -fastx -blast 1 -num_alignments 1 -v -no-best -other --paired_in -out2 -workdir "${FILE/_fwd_dedup.fq.gz/_sortmerna_p}"; echo "${FILE/_1.fastq.gz/} sortmerna end"; date; done

conda deactivate


# repair.sh on SortMeRNA's output
source activate bbmap_conda

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1.fastq.gz; do ACC="${FILE/_1.fastq.gz/}"; echo "${ACC} repair.sh start"; date; FWD="${ACC}_sortmerna_p/out/other_fwd.fq.gz"; REV="${ACC}_sortmerna_p/out/other_rev.fq.gz"; repair.sh in="$FWD" in2="$REV" -Xmx30g out="${ACC}_1_repaired_post_smrna.fastq.gz" out2="${ACC}_2_repaired_post_smrna.fastq.gz" overwrite=t; echo "${ACC} repair.sh end"; date; done

conda deactivate


# Trimmomatic on SortMeRNA's repaired output
source activate Trimmomatic

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1.fastq.gz; do ACC="${FILE/_1.fastq.gz/}"; echo "${ACC} Trimmomatic start"; date; FWD="${ACC}_1_repaired_post_smrna.fastq.gz"; REV="${ACC}_2_repaired_post_smrna.fastq.gz"; java -jar /home/paualico/.conda/envs/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE "$FWD" "$REV" "${ACC}_1_P_Trimmomatic.fastq.gz" "${ACC}_1_U_Trimmomatic.fastq.gz" "${ACC}_2_P_Trimmomatic.fastq.gz" "${ACC}_2_U_Trimmomatic.fastq.gz" -threads 25 ILLUMINACLIP:/home/paualico/.conda/envs/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25; echo "${ACC} Trimmomatic end"; date; done

conda deactivate

####
#### Kaiju ####
#### 30 cores

source activate kaiju

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1_P_Trimmomatic.fastq.gz; do FILE2=${FILE/_1_P_Trimmomatic.fastq.gz/_2_P_Trimmomatic.fastq.gz}; echo "${FILE/_1_P_Trimmomatic.fastq.gz/} kaiju start"; date; srun /home/paualico/.conda/envs/kaiju/bin/kaiju -t  /home/paualico/.conda/envs/kaiju/kaijudb/rvdb/nodes.dmp -f /home/paualico/.conda/envs/kaiju/kaijudb/rvdb/rvdb/kaiju_db_rvdb.fmi -i "$FILE" -j "$FILE2" -z 30 -E 0.00001 -v -o "${FILE/_1_P_Trimmomatic.fastq.gz/_kaiju_pair.txt}"; srun /home/paualico/.conda/envs/kaiju/bin/kaiju-addTaxonNames -t /home/paualico/.conda/envs/kaiju/kaijudb/rvdb/nodes.dmp -n /home/paualico/.conda/envs/kaiju/kaijudb/rvdb/names.dmp -r superkingdom,phylum,class,order,family,genus,species -u -i "${FILE/_1_P_Trimmomatic.fastq.gz/_kaiju_pair.txt}" -o "${FILE/_1_P_Trimmomatic.fastq.gz/_kaiju_pair_tax.txt}"; echo "${FILE/_1_P_Trimmomatic.fastq.gz/} kaiju end"; date; done

conda deactivate

####
#### ENSAMBLADO ####
####

module load singularity/3.4.1

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1.fastq.gz; do ACC="${FILE/_1.fastq.gz/}"; echo "${ACC} trinity start"; date; singularity exec -e /home/paualico/.conda/envs/trinity/trinityrnaseq.v2.14.0.simg Trinity --seqType fq --left "${ACC}_1_P_Trimmomatic.fastq.gz" --right "${ACC}_2_P_Trimmomatic.fastq.gz" --max_memory 100G --CPU 30 --full_cleanup --output "${ACC}_trinity_out_fp"; echo "${ACC} trinity end"; date; done

# E90N50 calculation
for FILE in /home/paualico/metavirome/FINAL/*/*/*Trinity.fasta; do python3 /home/paualico/samples/pau/script_E90N50.py --assembly "$FILE" --abundance "${FILE/_trinity_out_fp.Trinity.fasta/}_kallisto/abundance.tsv" --trinity_path /home/paualico/.conda/envs/trinity2/bin/ ;done > /home/paualico/metavirome/FINAL/E90N50stat.txt

####
#### ABUNDANCIA #### 
#### 20 cores

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1_P_Trimmomatic.fastq.gz; do ACC="${FILE/_1_P_Trimmomatic.fastq.gz/}"; FILE_2="${FILE/_1_P_Trimmomatic.fastq.gz/_2_P_Trimmomatic.fastq.gz}"; FILE_3="${FILE/_1_P_Trimmomatic.fastq.gz/_trinity_out_fp.Trinity.fasta}"; dir=${FILE/_1_P_Trimmomatic.fastq.gz/_kallisto}; echo "${ACC} kallisto start"; date; if [[ -f $FILE_2 ]]; then srun singularity exec -e /home/paualico/.conda/envs/trinity/trinityrnaseq.v2.13.2.simg /home/paualico/.conda/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/align_and_estimate_abundance.pl --seqType fq --left "$FILE" --right "$FILE_2" --transcripts "$FILE_3" --est_method kallisto --aln_method bowtie2 --thread_count 20 --trinity_mode --prep_reference --output_dir "$dir"; else srun singularity exec -e /home/paualico/.conda/envs/trinity/trinityrnaseq.v2.13.2.simg /home/paualico/.conda/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/align_and_estimate_abundance.pl --seqType fq --single "$FILE" --transcripts "$FILE_3" --est_method kallisto --aln_method bowtie2 --thread_count 20 --trinity_mode --prep_reference --output_dir "$dir"; fi; echo "${ACC} kallisto end"; date; done

####
#### DIAMOND ####
#### 30

source activate diamond

for FILE in /home/paualico/metavirome/FINAL/*/*/*_1.fastq.gz; do ACC="${FILE/_1.fastq.gz/}"; echo "${ACC} diamond start"; date; /home/paualico/.conda/envs/diamond/diamond blastx -d /home/paualico/.conda/envs/diamond/db/nr.dmnd -q "${ACC}_trinity_out_fp.Trinity.fasta" --very-sensitive -p 30 --index-chunks 1 -e 0.00001 --taxonlist 10239 -o "${ACC}_diamond_f.tsv" --outfmt 6 qseqid sseqid qstart qend sstart send evalue bitscore length pident staxids sscinames sskingdoms skingdoms sphylums; echo "${ACC} diamond end"; date; done

conda deactivate

####
#### SCRIPTS GENERADORES DE TABLAS ####
####

# Kaiju

python /home/paualico/metavirome/scripts/kaiju_compute_abundances.py /home/paualico/metavirome/FINAL/ _kaiju_pair_tax.txt

# Diamond

# get tpm per viral family for each accession
source activate taxonkit

for DMND in /home/paualico/metavirome/FINAL/*/*/*diamond_f.tsv; do KLST="${DMND/_diamond_f.tsv/_kallisto}"/abundance.tsv; python3 /home/paualico/metavirome/scripts/viral_contigs_abundance.py --diamond "$DMND" --kallisto "$KLST" --taxonkit /home/paualico/.conda/envs/taxonkit/bin/taxonkit --dmp_dir /home/paualico/.conda/envs/diamond/db/nr/ --verbose; done

conda deactivate

# add SRA ID to each composition
python /home/paualico/metavirome/scripts/add_sraID.py /home/paualico/metavirome/FINAL/

source activate R_env
# merge compositions into a single matrix
Rscript /home/paualico/metavirome/scripts/merging_tpm_tables.R /home/paualico/metavirome/FINAL/ /home/paualico/metavirome/FINAL/

conda deactivate



