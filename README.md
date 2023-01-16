# Varroa transcriptomics project 
Here I analyze transcriptomic sequence data of _Varroa destructor_, a parasite of the globally distributed western honey bee (_Apis mellifera_), which were collected from >50 regions of the world. Varroa mites were originally parasitizing the eastern honey bee (_Apis cerana_), but experienced a host switch after the global spread of _A. mellifera_ in the mid 20th century. Following this host switch, Varroa mites have successfully infested the western honey bee populations worldwide. 
This is a part of my PhD thesis project. 

## Files and directories 
### >data
Here you will find all the sample data including metadata. <br>
- "sample-list.xlsx" contains information on Varroa mite samples we have collected for the study. Includes collection site, date, host and parasite species. <br>
- "NCBI-phylo.xlsx" contains SRA accession from NCBI, virus name, collection year, region, and host species as well as samples from our study that was used for phylogenetics analysis. <br>
- "summary.xlsx" contains information of viruses used in variant calling (contig). NCBI SRA accession, length, name, abbreviation. <br>
- "viruses2020-ENA.fasta" is a fasta file containing all the virus seq used in variant calling as listed in data/viruses2020.xlsx

### > other files
- "cluster.json" is a file that dictates all the cluster configuration for each rule in snakemake <br>
- "Snakefile" is a file that dictates what rules and commands are run on snakemake <br>
