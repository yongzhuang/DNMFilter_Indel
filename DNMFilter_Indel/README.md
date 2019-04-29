# DNMFilter_Indel
# Introduction 
DNMFilter_Indel is a machine learning based tool designed to filter out false positive de novo indels obtained by any computational approaches from whole genome or exome sequencing data. It can be used as either a stand-alone tool to detect de novo indels or coupled with other commonly used de novo indel detection tool to improve specificity.
# Installation
The easiest way to get DNMFilter_Indel is to download the binary distribution from the DNMFilter_Indel github release page. Alternatively, you can build DNMFilter_Indel from source with gradle.
1. git clone --recursive https://github.com/yongzhuang/DNMFilter_Indel.git
2. cd DNMFilter_Indel/
3. gradle build  

If you want to run DNMFilter_Indel, you'll need:
1. Java SE Development Kit 8
2. R (Rscript exectuable must be on the path)
3. Runiversal package in R
4. gbm package (https://cran.r-project.org/web/packages/gbm/) in R  

You'll find the executable jar file in DNMFilter_Indel/build/libs/.
# File Instruction
1. bam list file (two columns, tab-separated)  
   Column 1: sampleID   
   Column 2: path of .bam file  

2. DNM file (.csv file, the first five columns are mandatory, comma-separated)  
   Column 1: familyID   
   Column 2: chromsome   
   Column 3: position  
   Column 4: reference allele  
   Column 5: alternative allele  
   Note: The DNM file should be first sorted by familyID and then by chromosome position.  

3. feature configuration file  
   The feature that values 1 is selected to train the model and filter DNMs, while the feature that values 0 is not selected.  

4. pedigree file (six columns, tab-separated)   
   Column 1: Family ID  
   Column 2: Individual ID  
   Column 3: Paternal ID  
   Column 4: Maternal ID  
   Column 5: Sex (1=male; 2=female; other=unknown)  
   Column 6: Phenotype  
   Note: The Family ID must be the same as the first column of DNM file, and the Individual ID must be the same as the first column of bam list file.
   
5. output file (.csv file, the six columns are mandatory, comma-separated)  
   Column 1: familyID   
   Column 2: chromsome   
   Column 3: position  
   Column 4: reference allele  
   Column 5: alternative allele
   Column 6: prediction score

# Contact 
   yongzhuang.liu@hit.edu.cn


