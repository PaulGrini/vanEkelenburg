# Step 3 generates the homozygous mapping filter. This filter is applied for Early stage marker (EE) and Late stage markers (ESR, DAL, TE1) 

- Make a folder */reference* that contains the polished Col-0 and Tsu-1 target sequences generated in step 1
- Prepare the directories for all replicates which contains softlinks to the trimmed reads using ```make_map_dirs.sh```
  Note that this is done separately for both sequencing dates
- Map the trimmed reads to the polished target sequences using ```grid.bowtie_map.sh``` which runs ```bowtie_map.sh``` in each subdirectory
- Extract and count informative read counts for each gene using ```grid.diff.sh``` which runs ```run_diff.sh``` in each subdirectory
- Obtain the genelists using ```run_filter_homozygous.sh```. This produces the homozygous mapping filter (SData 3) and generates the genelists used in step 7
