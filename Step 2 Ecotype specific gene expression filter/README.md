# Step 2 generates the ecotype specific gene expression filter. This filter is applied for Early stage marker (EE) and Late stage markers (ESR, DAL, TE1) 

- Make a folder */reference* that contains the polished Col-0 and Tsu-1 target sequences generated in step 1
- Prepare the directories for all replicates which contains softlinks to the trimmed reads using ```make_map_dir.sh```
  Note that this is done separately for both sequencing dates
- Map the trimmed reads to the polished target sequences using ```grid.bowtie_map.sh``` which runs ```bowtie_map.sh``` in each subdirectory
- Extract read counts for each gene using ```grid.samstats.sh``` which runs ```samstats.sh``` in each subdirectory
```diff
! Important, rename the counts.tsv files manually so that they contain the replicate number,
! ecotype and sequencing date: 01-Col-0-Early-BR1.20190705.tsv, make a directory called */Counts* 
! and collect all tsv files in this directory 
```
- Make another directory called '/Merged_read_counts'
- Prepare the read files for differential gene expression with ```Preparing_homozygous_reads_Early.R``` and ```Preparing_homozygous_reads_Late.R```
> These scripts do a number of things to prepare the files used for differential gene expression:
With bowtie, reads were mapped to either the Col-0 or Tsu-1 allele. Here we don't need this and therefore read counts of the same gene are summated.
Read counts from both sequencing dates for the same sample are combined

- Make a directory called */Output* and run DESeq2 using ```DESeq2_homozygous_reads_Early.R``` and ```DESeq2_homozygous_reads_Late.R```
> This produces MAplots (SFigure 4)

- In order to obtain the genelists for the ecotype specific gene filter, run ```Generating_ecotype_specific_filter.R```
> This merges DESeq2 input and output (SData 4) and generates the genelists used in step 7
