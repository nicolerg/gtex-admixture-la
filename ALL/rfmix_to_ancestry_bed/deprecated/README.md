# GTEx v8 Local Ancestry Calls  
## Contact:  
Nicole Gay  
Montgomery Lab   
nicolerg@stanford.edu    

# Files  
This package contains two files:  
1. **rfmix.viterbi.master.tsv.gz**  
2. **master.ancestry.bed**  

## rfmix.viterbi.master.tsv.gz  
This file contains SNP-wise ancestry assignments as determined by RFMix v1.5.4. Variant information is given in the first four columns. The remaining columns provide the corresponding ancestry assigments for each haplotype of 114 admixed individuals, defined as self-reported African American and Asian subjects in GTEx v8. The two haplotypes for each individual are labelled `SUBJECT-ID_A` and `SUBJECT-ID_B`. The possible ancestry assignments are 1, 2, and 3. These correspond to Asian, European, and African populations, respectively. More specifically, the 1000 Genomes reference panels used to define each population are as follows:  
1. Asian: CHB, JPT  
2. White: CEU  
3. African: YRI, GWD, MSL, ESN  
Variants for this analysis were limited to those present in both GTEx v8 and 1000 Genomes. Coordinates are in hg19 for compatibility with the required input format for RFMix v1.5.4 (variant positions must be given in cM, not bp, and there is no genetic map for hg38).   

## master.ancestry.bed  
This BED file provides a summarized version of the information in **rfmix.viterbi.master.tsv.gz**. Given the ancestry assignments in **rfmix.viterbi.master.tsv.gz** and the posterior probabilities of each population at each SNP (another output of RFMix, not given here), consecutive variants with the same local ancestry assignment (posterior probability > 0.9) are collapsed into blocks on a per-haplotype basis. Note that there is always a small range of bases with unknown ancestry between blocks with different local ancestry. Rows are ordered first by haplotype, then by chromosome. Haplotypes from the same individual follow each other in the row ordering. Subject-level BED files were generated using a script written by Alicia Martin (https://github.com/armartin/ancestry_pipeline/blob/master/collapse_ancestry.py).  
