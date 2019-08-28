# genarch_malechoice_Drosophila
Repository of scripts used to align whole genome sequence and estimate ancestry for DOI: 

**Sequence analysis pipeline in order**

1. Align all files to reference (paired-end reads) and generate .bam files:

- BWA_PEaln_SAMPLE.sh

2. Call SNPs for parent strains: 

- SNPs_parents2L.sh  
- SNPs_parents2R.sh  
- SNPs_parents3L.sh  
- SNPs_parents3R.sh
- SNPs_parentsX.sh     
- SNPs_parents4.sh  

3. Compile list of parent SNPs from .vcf file (custom python scripts):

- SNPs_parents2L.py
- SNPs_parents2R.py
- SNPs_parents3L.py
- SNPs_parents3R.py
- SNPs_parentsX.py
- SNPs_parents4.py  

4. Generate pileup files for each chromosome arm of individuals from mapping population:

- SAMPLE_pileup2L.sh
- SAMPLE_pileup2R.sh
- SAMPLE_pileup3L.sh
- SAMPLE_pileup3R.sh
- SAMPLE_pileupX.sh
- SAMPLE_pileup4.sh

5. Sort pileups for sliding window analysis (custom python scripts):

- SAMPLE_sortpileup2L.py
- SAMPLE_sortpileup2R.py
- SAMPLE_sortpileup3L.py
- SAMPLE_sortpileup3R.py
- SAMPLE_sortpileupX.py
- SAMPLE_sortpileup4.py

6. Sliding window ancestry:

- Ancestry2L_SAMPLE.py
- Ancestry2R_SAMPLE.py
- Ancestry3L_SAMPLE.py
- Ancestry3R_SAMPLE.py
- AncestryX_SAMPLE.py
- Ancestry4_SAMPLE.py

**NOTE: For ‘SAMPLE’ scripts, we used a custom python script to iteratively generate files that replaced ‘SAMPLE’ with each of the mapping population’s sample numbers. Example file provided:

- AncestryX.py
