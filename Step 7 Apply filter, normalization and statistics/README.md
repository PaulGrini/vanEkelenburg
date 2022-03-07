# Step 7 applies filters, normalizations and statistics on the IRP output (Step 4)

- Collect all generated .tsv files from step 4 to this location and make sure that the normalization factors are included in **ColTsu.model.tsv**
- Collect the generated filter list files (step 2 and 3) 
> Early_ecotype_pass_filter.tsv and Late_ecotype_pass_filter.tsv (Step 2)
 Early_genes_pass_filter.tsv and Late_genes_pass_filter.tsv (Step 3)

- Perform filtering, normalization and statistics on informative reads from markers using ```run_stats.sh```

> The ouput will be used in Step 8 Data analysis
