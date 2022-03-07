# Step 5 extracts library normalization factors for each replicate of the marker lines

- Make a folder */reference* that contains the polished Col-0 and Tsu-1 target sequences generated in step 1
- Prepare the directories for all replicates which contains softlinks to the trimmed reads using ```make_map_dirs.sh```
- Sort the Primary.bam files using ```SortedBAM.sh```
- Count the mapped reads using ```grid.samstats.sh``` which runs ```samstats.sh``` in each subdirectory 
- Rename the counts.tsv files to their respective sample/replicate name using ```Rename_counts.sh```

> Proceed to the directory 'Normalization factors'
- Copy all renamed counts.tsv files using ```copy.sh``` and remove unnecessary columns using ```Remove_column.sh```
- Prepare the files for the next step using ```prepare_heterozygous.sh``` which collates individual replicates in one file names as follows: 'pos.EE1.counts.csv' (for EE)
- Extract Normalization factors using ```Normalization_factor.sh``` that are found in files called 'EE1pos.normfactors.csv' (for EE)
> Normalization factors need to be manually written inside ColTsu.model.tsv located in 'Step 7 Apply filter, normalization and statistics' under Norm
