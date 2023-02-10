# *PhyloRooting*

Welcome to PhyloRooting repository! Here you will find the scripts to infer phylogenomic rooting neighborhood as described in blablah. 

The main analysis is done with pgroot.py. The input of the analysis is a concatenated file of rooted trees whose branch values represent AD (run MAD with -u flag)
The program will output two files:

`*.df` contains a dataframe with the AD values per tree and partition
`*.mat` describes the OTUs composition of the candidate root partitions. 

If the flag `--neighborhood` is used, the inference of a root neighborhood will be performed. 

### Example datasets

The data used in the study blah can be additionally found here. The analysis on the proteobacteria dataset would run as:
`pgroot.py -t Proteobacteria.nwk.AD_unrooted -AD Proteobacteria.df -p Proteobacteria_part.mat --neighborhood`
