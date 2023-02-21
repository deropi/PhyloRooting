# *PhyloRooting*

Welcome to PhyloRooting repository! Here you will find the scripts to infer phylogenomic rooting neighborhood as described in https://doi.org/10.1101/758581. 

The main analysis is done with pgroot.py. The input of the analysis is a concatenated file of rooted trees whose branch values represent AD (run MAD with -u flag). 

Sequence headers ***must follow the following format:***

`>12_1`

`MDVGKKKTKGC`

`>12_2`

`MDVGKKKTKGC`

`>14_1`

`MDVGKKKTKGC`

Where the first element corresponds to the genome/species id and the second element corresponds to the copy number. For instance, in this example, species 12 contains two paralogs for the same protein, while species 14 has only one. 

The program will output two files:

`*.df` contains a dataframe with the AD values per tree and partition
`*.mat` describes the OTUs composition of the candidate root partitions. 

If the flag `--neighborhood` is used, the inference of a root neighborhood will be performed. 

### Example datasets

The data used in the study blah can be additionally found here. The analysis on the proteobacteria dataset would run as:
`pgroot.py -t Proteobacteria.nwk.AD_unrooted -AD Proteobacteria.df -p Proteobacteria_part.mat --neighborhood`
