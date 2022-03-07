Polishing of Col-0 and Tsu-1 reference transcriptomes using pilon.

Before running pilon, make the subdirectories for each of the reads by running ```make_map_dir.sh```, fill in the location of the trimmed reads.

Continue with making consensus_0, which originates from Araport11 using ```Initial_mapping_to_Araport11.sh```
The araport11 reference that was used is located in **reference**, where it has been renamed to Araport.fasta..

Proceed with mapping to consensus_0, by running 
```
grid.bowtie_map.sh 0 # The 0 specifies round 0
```
  This will generate directories named **map_to_consensus_0** and this contains subdirectories of all reads that mapped

Follow the first mapping round with 
```
grid.pilon.sh 1 # The 1 specifies the generation of consensus reference 1.
```
  The number of changes made to a reference transcriptome compared to the previous is highlighted by a *.changes* file.

Repeat the procedure of mapping to consensus 1 by running 
```
grid.bowtie_map.sh 1
grid.pilon.sh 2 # This generates consensus reference 2
# etc
```

Continue with the process of mapping and polishing until the number of changes has dropped to a desired level. 

After polishing with pilon is completed, merge the polished reference transcriptomes by running ```merge_final_consensus_references.sh```
  Include the path to the consensus_ directory that was last generated
