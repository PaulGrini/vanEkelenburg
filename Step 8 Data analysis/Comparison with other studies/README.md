# Comparison of identified imprinted genes with previous studies

- Other studies only have gene codes, and not transcript number. Therefore transcriptnumbers first need to be removed.
> Copy informative reads with gene description as prepared in */Data analysis* using ```Prepare_heterozygous_reads.R``` 
Use the terminal to remove transcript numbers with
```
sed 's/\(AT[0-9]G[0-9]*\)\.[0-9]*/\1/g' EE_inf_reads_and_gene_description.csv > EEpos_no_transcript.csv
sed 's/\(AT[0-9]G[0-9]*\)\.[0-9]*/\1/g' ESR_inf_reads_and_gene_description.csv > ESRpos_no_transcript.csv
sed 's/\(AT[0-9]G[0-9]*\)\.[0-9]*/\1/g' TE1_inf_reads_and_gene_description.csv > TE1pos_no_transcript.csv
```

- Compare the overlap of identified imprinted genes using ```Comparison_with_other_studies.R```
> Genes that were identified as MEGs/PEGs are provided in 'Other_studies2.tsv'
  Genes that were identified in this study and also in a previous study were extracted (SData 8)
  This script also prepares the input to generate the barplot showing the increase in overlap after removing genes not analyzed in this study (SData 9)
  The overlap before and after excluding the genes that were not analyzed in this study is provided in SData 9

- Generate the barplot using ```Overlap_barplot.R```
> 'SData 9 Overlap in imprinted genes when only looking at genes present in our data.xlsx' needed to have an extra sheet for Rstudio with the following structure:
  Column A: Domains (EE, ESR, TE1, All); Column B: Study (Picard, Hornslien, Pignatta, Del_Toro); Column C: Imprinting (MEGs, PEGs); 
  Column D: Overlap (link to the overlap percentages on each sheet); Column E: Type (All, Filtered) 
  This was named 'SData 9 Overlap in imprinted genes when only looking at genes present in our data2.xlsx' and is located in */Input*
