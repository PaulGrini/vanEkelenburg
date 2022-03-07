```Trimgalore.sh``` is used to run the program that will trim the read.
```Grid.trimgalore.sh``` is a shell script that calls on ```trimgalore.sh``` for every read inside the directory.

```Trimgalorehetero.sh``` and ```grid.trimgalorehetero.sh``` are used to trim reads originating from the marker lines. 
These reads have been sequenced with the ISPCR primer and illumina adapters, so they both need to be trimmed. Adapters and primers are placed in the file named ```adapters.fa```. 
