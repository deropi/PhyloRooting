# *PhyloRooting*

Welcome to PhyloRooting repository! Here you will find the scripts to infer phylogenomic rooting neighborhood as described in https://doi.org/10.1101/758581. 

The main analysis is done with pgroot.py. The input of the analysis is a concatenated file of rooted trees whose branch values represent AD (run MAD with -u flag). 

##Important! Sequence headers ***must follow the following format:***

`>12_1`

`MDVGKKKTKGC`

`>12_2`

`MDVGKKKTKGC`

`>14_1`

`MDVGKKKTKGC`

Where the first element corresponds to the genome/species id and the second element after the underscore indicates the copy number. For instance, in this example, species 12 contains two paralogs of the same protein, while species 14 has only one. 

The program will output two files:

`*.df` contains a dataframe with the AD values per tree and partition
`*.mat` describes the OTUs composition of the candidate root partitions. 

If the flag `--neighborhood` is used, the inference of a root neighborhood will be performed. 

### Example datasets

Additionally, the data sets used in the study can be found here. The analysis for the proteobacteria dataset can be run as follows:

`pgroot.py -t Proteobacteria.nwk.AD_unrooted -AD Proteobacteria.df -p Proteobacteria_part.mat --neighborhood`
