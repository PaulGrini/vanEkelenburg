# Step 8 performs all the data analysis

- Analysis of qPCR data is performed in */qPCR*

- Copy the **.filtered.final.csv** files generated in step 7 to the 'Input' directory
- Make the directory */Total inf read and ecotype specific gene filter*
- Make the directory */Output*

- Distribution of genes (imprinted, parental bias, non-parental bias, not enough informative reads) is performed in */Distribution of genes*
>  This generates the piecharts (SFigure 5)

- Prepare the generated files for data analysis using ```Prepare_heterozygous_reads.R```
> This script renames columns, removes unnecessary columns, includes a column containing the sum of all informative reads and adds gene descriptions for each gene
  The output is saved in 'Output' (SData 5)
  Then the script filters on >30 total informative reads to ensure enough depth of informative reads
- Imprinting analysis for each marker line is performed using ```Informative_read_differential_expression_Col_vs_Tsu.R```
> MEGs and PEGs are selected based on positive or negative log2FC respectively and adjusted p-value < 0.05
  All identified MEGs and PEGs are placed in one excel file (SData 6)
  Furthermore, MEGs and PEGs unique for each marker and overlapping are identified (SData 10; Figure 4)
  
- The overlap between MEGs and PEGs identified in this study and previously published datasets is analyzed in */Comparison with other studies*

- Spatial and temporal specific imprinting is analyzed in */Log2 fold change plots of informative reads*
> Temporal specific imprinting:
  Fold changes of imprinted MEGs/PEGs in EE, TE1 and both EE and TE1 were collected (files made by ```Prepare_heterozygous_reads.R```)
  Only genes that had sufficient informative reads in both markers were analyzed and the other genes were excluded (SData 11)
  Fold changes were plotted against each other (Figure 5) and input for this plot was extracted (SData 11)
  
> Spatial specific imprinting:
  Fold changes of imprinted MEGs/PEGs in ESR, TE1 and both ESR and TE1 were collected (files made by ```Prepare_heterozygous_reads.R```)
  Only genes that had sufficient informative reads in both markers were analyzed and the other genes were excluded (SData 12)
  Fold changes were plotted against each other (Figure 6) and input for this plot was extracted (SData 12)
