1) virfinder-parallel-github_v3_25cores.R runs VirFinder in parallel for the obtention of p-values (H_0: non viral, H_1: viral) for each input sequence. Modification of https://github.com/rec3141/VirFinder/blob/master/linux/VirFinder/R/parVF.pred.R

2) postproc_vf.R adds a new column for FDR adjustment of all p-values. Another column is formed by a logical vector TRUE when FDR < 0.05 and FALSE otherwise. A third logical vector column is added to identify observations in the percentile 30 lowest p-values.

3) virfinder2fasta.py parses the postprocessed file and selects contigs with  FDR < 0.05 or 30% lowest p-values (separate outputs) from the complete fasta file, yielding the final subset of contigs classified as viral by VirFinder.



