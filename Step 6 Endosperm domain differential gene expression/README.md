# Step 6 performs differential gene expression analysis between marker lines

- Make a folder */reference* that contains the polished Col-0 and Tsu-1 target sequences generated in step 1
- Prepare the directories for all replicates which contains softlinks to the trimmed reads using ```make_map_dir.sh```
  Note that this is done separately for both sequencing dates
- Map the trimmed reads to the polished target sequences using ```grid.bowtie_map.sh``` which runs ```bowtie_map.sh``` in each subdirectory
- Extract read counts for each gene using ```grid.samstats.sh``` which runs ```samstats.sh``` in each subdirectory
- Rename counts.tsv files using ```rename_counts.sh``` and copy the tsv files in 'input' using ```copy.sh```
> Make a directory */Genefilter* and place the Early_genes_pass_filter.tsv and Late_genes_pass_filter.tsv files generated in step 3
  This excludes the genes that did not pass the homozygous gene filter list in the next step
- Prepare the files for differential gene expression analysis using ```Prepare_heterozygous_reads.R```
- Perform differential expression with DESeq2 using ```DESeq2_heterozygous_reads.R```
> This generates the PCA plot (Figure 2A) from which was determined to exclude DAL from further analyses
DESeq2 is repeated without DAL to produce MAplots (Figure 2A, SFigure 2A & 2B)
Make a directory 'Output' to collect the DESeq2 output data (SData 2)
