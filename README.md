# Varroa transcriptomics project 
Here we analyze transcriptomic sequence data of _Varroa destructor_, a parasite of the globally distributed western honey bee (_Apis mellifera_), which were collected from >50 regions of the world. Varroa mites were originally parasitizing the eastern honey bee (_Apis cerana_), but experienced a host switch after the global spread of _A. mellifera_ in the mid 20th century. Following this host switch, Varroa mites have successfully infested the western honey bee populations worldwide. <br>
Using transcriptomic data, viruses affiliated with the host honey bee and/or the parasite mite were mined from over 200 mite samples, collected from 50+ regions of the world, over 4 decades between 1980s-2010s. In this [PNAS paper]([https://doi.org/10.1073/pnas.2301258120]), we examined varroa and honey bee associated viruses, and mapped out its distribution. We also performed phylogenetic analyses which reveal that DWV-a originated in Asia. 
Raw sequence data are available under DDBJ/NCBI BioProject [[PRJDB14940](https://ddbj.nig.ac.jp/resource/bioproject/PRJDB14940)]. <br>
This project is a part of my PhD thesis at Okinawa Institute of Science and Technology. <br>

## Files and directories 
### data
Here you will find all the sample data including metadata. <br>
- "sample-list.xlsx" contains information on Varroa mite samples we have collected for the study. Includes collection site, date, host and parasite species. <br>
- "NCBI-phylo.xlsx" contains SRA accession from NCBI, virus name, collection year, region, and host species as well as samples from our study that was used for phylogenetics analysis. <br>
- "summary.xlsx" contains information of viruses used in variant calling (contig). NCBI SRA accession, length, name, abbreviation. <br>
- "viruses2020-ENA.fasta" is a fasta file containing all the virus seq used in variant calling as listed in data/viruses2020.xlsx
- "DWV-A.xml" and "DWV-B.xml" contain virus alignments used in the phylogeography analysis

### other files
- "cluster.json" is a file that dictates all the cluster configuration for each rule in snakemake <br>
- "Snakefile" is a file that dictates what rules and commands are run on snakemake <br>

## Virome analysis
- In the virome folder, there is an r script and data files used to generate the relevant figures.

## Contact
If you have any questions, please reach Nonno Hasegawa, at <br>
nonno.hasegawa@outlook.com

## Preprint
The preprint is available on [biorxiv](https://www.biorxiv.org/content/10.1101/2023.01.21.525007v2.abstract).

## Manuscript accepted and published [here](https://doi.org/10.1073/pnas.2301258120).
