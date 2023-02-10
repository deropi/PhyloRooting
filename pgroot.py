#!/usr/bin/env python3
#----------------------- LOADING MODULES -----------------------------------------------

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString, _count_clades, _tree_to_bitstrs, _clade_to_bitstr
import sys
import pandas as pd
import argparse
import itertools
import scipy.stats as stats
import numpy as np
from datetime import datetime

#---------------------- PARSING ARGUMENTS ---------------------------------------------

parser = argparse.ArgumentParser() # Generate the parser

parser.add_argument('-t',metavar='trees',help='This file contains concantenated unrooted Newick trees where AD values are shown as branch label (pgroot output from MAD)')

parser.add_argument('-AD',metavar='AD_df',help='CSV file containing a dataframe with AD values per partition and tree')

parser.add_argument('-p',metavar='partition_file',help='CSV file with partitions information')

parser.add_argument('--neighborhood', action='store_true')

args = parser.parse_args()

#--------------------- FUNCTIONS -----------------------------------------------------
def _tree_to_bitstrs(tree):
	"""Create a dict of a tree's clades to corresponding BitStrings (PRIVATE)."""
	clades_bitstrs = {}
	term_names = [term.name for term in tree.find_clades(terminal=True)]
	for clade in itertools.chain(tree.find_clades(terminal=True), tree.find_clades(terminal=False)):
		bitstr = _clade_to_bitstr(clade, term_names)
		clades_bitstrs[clade] = bitstr
	return clades_bitstrs

def get_bitstrings(tree):
	terminals = [x.name.split('_')[0] for x in tree.get_terminals()]
	clade_bitstrings = _tree_to_bitstrs(tree)
	return (terminals, clade_bitstrings)

def append_partition(diction, clade, count, part_dict, bitstring, genfam_type, AD='Non-Inf'):
	if AD != None:
		diction['Tree'].append(count)
		diction['AD'].append(AD)
		diction['Partition_id'].append(part_dict[bitstring])
		diction['Partition'].append(bitstring)
		diction['GeneFam'].append(genfam_type)

def add_AD(AD,diction,tree):
	if AD != None:
		diction[tree].append(float(AD))

#--------------------- PROCESS --------------------------------------------------------
start_time = datetime.now()
# Reading and classifying trees into CSC, CMC, PSC and PMC 
otus = []
treeclas = {}
trees = list(Phylo.parse(args.t, "newick"))
print ('Parsing trees ...\n')
count = 0
for tree in trees:
	count += 1
	termi = tree.get_terminals()
	names = [x.name.split('_')[0] for x in termi]
	treeclas[count] = [names, tree, str()]
	otus +=names
otus_set = set(otus)
for tree in treeclas:
	if set(treeclas[tree][0]) == otus_set and len(treeclas[tree][0]) == len(otus_set):
		treeclas[tree][2] = 'CSC'
	elif set(treeclas[tree][0]) == otus_set and len(treeclas[tree][0]) > len(otus_set):
		treeclas[tree][2] = 'CMC'
	elif set(treeclas[tree][0]) != otus_set and len(treeclas[tree][0]) > len(set(treeclas[tree][0])):
		treeclas[tree][2] = 'PMC'
	elif set(treeclas[tree][0]) != otus_set and len(treeclas[tree][0]) == len(set(treeclas[tree][0])):
		treeclas[tree][2] = 'PSC'
CSC,CMC,PSC,PMC = [x for x in treeclas if treeclas[x][2]=='CSC'], [x for x in treeclas if treeclas[x][2]=='CMC'], [x for x in treeclas if treeclas[x][2]=='PSC'], [x for x in treeclas if treeclas[x][2]=='PMC']
if len(CSC)==0:
        raise ValueError("No Complete Single-copy tree has ben found.")
print ('Parsed:\n' + str(len(otus_set)) + ' species\n' + str(count) + ' trees\n\n\tSingle-copy\tMulti-copy\nComplete\t'+ str(len(CSC)) + '\t' + str(len(CMC)) + '\nPartial\t' + str(len(PSC)) + '\t' + str(len(PMC)))

## CSC trees analysis and candidate root partitions identification

# Terminal names are set to the order in the first CSC provided tree
name_order = trees[CSC[0]].get_terminals()
names = [x.name.split('_')[0] for x in name_order]
if not [len(x.name.split('_')) for x in name_order].count(2) == len(names):
	raise ValueError("OTU names must have two components: X_1. Where X is the species id and the number after '_' represents the copy number.")
ADs_pertree = {}
diction_csc = {'Tree':[], 'Partition':[], 'Partition_id':[], 'AD':[], 'GeneFam':[]}
part_dict = {}
for tree_id in CSC:
	tree = treeclas[tree_id][1] 
	ADs_pertree[tree_id] = []
	names_new, clade_bitstrs = get_bitstrings(tree)
	index = [names_new.index(i) for i in names] #position in names_new which corresponds to original order
	for clade in clade_bitstrs:
		AD = clade.branch_length
		add_AD(AD,ADs_pertree,tree_id)
		bitstr_ori = clade_bitstrs[clade]
		bitstr = ''.join([bitstr_ori[i] for i in index])
		bitstrcomp = ''.join('1' if x == '0' else '0' for x in bitstr)
		if bitstr in part_dict:
			bitstring = bitstr
		elif bitstrcomp in part_dict:
			bitstring = bitstrcomp
		else:
			part_dict[bitstr] = len(part_dict.values())+1
			bitstring = bitstr
		if len(set(list(bitstring))) > 1:
			append_partition(diction_csc, clade, tree_id, part_dict, bitstring, 'CSC', AD)

# Storing tree_data into AD_df
csc_part = pd.DataFrame(data = diction_csc)
csc_root_candidates = csc_part[csc_part.groupby("Tree")["AD"].transform('min') == csc_part['AD']].copy()
csc_root_candidates['Counts'] = csc_root_candidates.groupby('Tree')['Partition_id'].transform('count')
csc_root_candidates['consensus'] = 1/csc_root_candidates.Counts
csc_root_candidates['Voted_root'] = csc_root_candidates.groupby('Partition_id')['consensus'].transform('sum')
csc_root_candidates['MAD_root'] = 100*(csc_root_candidates.Voted_root/len(csc_root_candidates.Tree.unique()))
a = csc_root_candidates[['Partition_id', 'Voted_root', 'MAD_root']].drop_duplicates().sort_values(by='MAD_root', ascending=False)
a = a.reset_index()
a = a.rename(columns={"index":"Part_id"})
a['Part_id'] = a.index + 1
print ('\nCandidate root partition counts among Complete Single-Copy trees:\n')
print (a[['Part_id', 'Voted_root', 'MAD_root']])
print ('\nVoted_root: number of CSC trees where this was the best root partition\nMAD_root: percentage of CSC trees where this was the best root partition\n')
a = a.drop(['Voted_root', 'MAD_root'], axis=1)
b  = csc_root_candidates.merge(a, on=['Partition_id'], how = 'left').drop_duplicates(subset=['Tree', 'Partition_id'])
b = b.drop(['Partition_id', 'Voted_root', 'MAD_root', 'Counts', 'consensus'],axis=1)
b = b.rename(columns={'Part_id': 'Partition_id'})
part_dict = dict(zip(b.Partition, b.Partition_id))
part_dict_rev = dict(zip(b.Partition_id, b.Partition))
diction = b.to_dict('list')
root_candidates = sorted(part_dict_rev.keys())
print ('\n'+str(len(part_dict)) + ' partitions are identified.')

# Storing partitions information into p
with open(args.p, 'w') as out:
	header = ','.join(['Partition_id'] + names )
	out.write(header + '\n')
	for part in root_candidates:
		row = str(part) + ','+ ','.join(list(part_dict_rev[part]))
		out.write(row + '\n')

## CMC trees analysis

print ('\nIdentifying candidate root partitions in non-CSC trees...')
for tree_id in CMC:
	tree = treeclas[tree_id][1]
	ADs_pertree[tree_id] = []
	names_new, clade_bitstrs = get_bitstrings(tree)
	for clade in clade_bitstrs:
		AD =clade.branch_length
		add_AD(AD,ADs_pertree,tree_id)
		bitstr = clade_bitstrs[clade]
		dupl = {}
		for i in range(0, len(names_new)):
			if not names_new[i] in dupl:
				dupl[names_new[i]] = list()
			dupl[names_new[i]].append(bitstr[i])
		dupl_splits = [len(set(x)) for x in dupl.values()]
		if len(set(dupl_splits)) == 1:
			names_new_clean = []
			bitstr_clean = []
			for e, f in zip(names_new, bitstr):
				if not e in names_new_clean:
					names_new_clean.append(e)
					bitstr_clean.append(f)
			index = [names_new_clean.index(i) for i in names]
			bitstr2 = ''.join([bitstr_clean[i] for i in index])
			bitstrcomp = ''.join('1' if x == '0' else '0' for x in bitstr2)
			if bitstr2 in part_dict:
				append_partition(diction, clade, tree_id, part_dict, bitstr2, 'CMC', AD)
			elif bitstrcomp in part_dict:
				append_partition(diction, clade, tree_id, part_dict, bitstrcomp, 'CMC', AD)
print ('\nCMC complete...')

## PSC trees analysis

for tree_id in PSC:
	tree = treeclas[tree_id][1]
	ADs_pertree[tree_id] = []
	names_new, clade_bitstrs = get_bitstrings(tree)
	index = [names_new.index(i) for i in names if i in names_new]
	names_new_sorted = [names_new[i] for i in index]
	for clade in clade_bitstrs:
		AD =clade.branch_length
		add_AD(AD,ADs_pertree,tree_id)
		bitstr = clade_bitstrs[clade]
		bitstr2 = ''.join([bitstr[i] for i in index])
		bitstrcomp = ''.join('1' if x == '0' else '0' for x in bitstr2)
		for each in part_dict:
			index_reduced = [names.index(i) for i in names_new_sorted]
			reduced_bits = ''.join([each[i] for i in index_reduced])
			if (bitstr2 == reduced_bits or bitstrcomp == reduced_bits) and len(set(list(reduced_bits))) >1:
				append_partition(diction, clade, tree_id, part_dict, each, 'PSC', AD)
			elif (bitstr2 == reduced_bits or bitstrcomp == reduced_bits) and len(set(list(reduced_bits))) == 1:
				append_partition(diction, clade, tree_id, part_dict, each, 'PSC')
print ('\nPSC complete...')

## PMC trees analysis

for tree_id in PMC:
	tree = treeclas[tree_id][1]
	ADs_pertree[tree_id] = []
	names_new, clade_bitstrs = get_bitstrings(tree)
	for clade in clade_bitstrs:
		AD =clade.branch_length
		add_AD(AD,ADs_pertree,tree_id)
		bitstr = clade_bitstrs[clade]
		dupl = {}
		for i in range(0, len(names_new)):
			if not names_new[i] in dupl:
				dupl[names_new[i]] = list()
			dupl[names_new[i]].append(bitstr[i])
		dupl_splits = [len(set(x)) for x in dupl.values()]
		if len(set(dupl_splits)) == 1:
                #then, remove duplicate
			names_new_clean = []
			bitstr_clean = []
			for e, f in zip(names_new, bitstr):
				if not e in names_new_clean:
					names_new_clean.append(e)
					bitstr_clean.append(f)
			index = [names_new_clean.index(i) for i in names if i in names_new_clean]
			names_new_clean_sorted = [names_new_clean[i] for i in index]
			bitstr2 = ''.join([bitstr_clean[i] for i in index])
			bitstrcomp = ''.join('1' if x == '0' else '0' for x in bitstr2)
			for each in part_dict:
				index_reduced = [names.index(i) for i in names_new_clean_sorted]
				reduced_bits = ''.join([each[i] for i in index_reduced])
				if (bitstr2 == reduced_bits or bitstrcomp == reduced_bits) and len(set(list(reduced_bits))) >1:
					append_partition(diction, clade, tree_id, part_dict, each, 'PMC', AD)
				elif (bitstr2 == reduced_bits or bitstrcomp == reduced_bits) and len(set(list(reduced_bits))) == 1:
					append_partition(diction, clade, tree_id, part_dict, each, 'PMC')
print ('\nPMC complete...')

## Filling informative tree partitions with maxAD

ADmax = {key:max(value) for (key, value) in ADs_pertree.items()}
df = pd.DataFrame(data = diction)
df_unstack = pd.DataFrame(df.groupby(['Tree','Partition_id']).AD.first().unstack())
df_unstack = df_unstack.reset_index()
df_melt = pd.melt(df_unstack, id_vars='Tree', value_vars=list(df_unstack.columns[1:]), var_name= 'Partition_id', value_name='AD')
df_melt.AD = df_melt.AD.fillna(df_melt.Tree.map(ADmax))
df_melt.AD = pd.to_numeric(df_melt.AD, errors='coerce') 
df_ready = df_melt.merge(df[['GeneFam', 'Tree']], on=['Tree'], how = 'left').drop_duplicates(subset=['Tree', 'Partition_id'])
df_ready.to_csv(args.AD, index=False, na_rep='nan')
print ('\nADs and partitions stored in ' + args.AD + ' and ' + args.p )

## Root neighborhood

#df.to_csv('nonmaxAD_' +  args.AD, index=False, na_rep='nan')
if args.neighborhood == True:
	print ('\nProceeding with root neighborhood inference...\n')
else:
	print ('\n**COMPLETED**\nRuntime: ' + str(datetime.now() - start_time))
	exit()
# Convert df back into dictionaries

df_ready = df_ready.dropna()
dict_all = {k: f.groupby('Tree')['AD'].apply(float).to_dict() for k, f in df_ready.groupby('Partition_id')}
df_worst = df[(df['AD'] != 'Non-Inf')]
dict_worst = {k: f.groupby('Partition_id')['AD'].apply(float).to_dict() for k, f in df_worst.groupby('Tree')}

initial = len(dict_all)
iter = 1
tests = 0
print ('Iteration\tRejected partition\tFDR-pvalue\t#trees')
while len(dict_all)>1 and (initial > len(dict_all) or iter == 1):
	initial = len(dict_all)
	pvalues = {}
	for part in list(dict_all):
		part_AD = list()
		worst_AD = list()
		for tree in dict_all[part]:
			if tree in dict_worst:
				subdict = dict_worst[tree]
				values = [subdict[p] for p in subdict if p != part]
				if len(values) > 0:
					part_AD.append(dict_all[part][tree])
					worst_AD.append(max(values))
		substracted = [x-y for x,y in zip(part_AD,worst_AD)]
		if len(substracted) >=11:
			wilcox = stats.wilcoxon(substracted,alternative='greater', nan_policy='omit')
			pvalues[part] = [wilcox.pvalue, '', len(substracted)]
	pvals = [x[0] for x in pvalues.values()]
	pvals.sort()
	ranks = range(1, len(pvals)+1)
	tests += len(pvals)
	for i in range(0,len(pvals)):
		p = pvals[i]
		rank = ranks[i]
		fdr_pval = p*tests/rank
		for key in pvalues:
			if pvalues[key][0] == p and pvalues[key][1] == '':
				pvalues[key][1] = fdr_pval
	fdr_sorted = [x[1] for x in pvalues.values()]
	mini = min(fdr_sorted)
	rej = ''
	for key in pvalues:
		if pvalues[key][1] == mini:
			rej = key
	print (str(iter) + '\t' + str(rej) + '\t' + str(pvalues[rej][1]) + '\t' + str(pvalues[rej][2]))
	if mini <= 0.05:
		del dict_all[rej]
		for tree in dict_worst:
			if rej in dict_worst[tree]:
				del dict_worst[tree][rej]
	iter += 1

print ('\nRetained partition(s): ' + str(list(dict_all.keys())))
print ('\n**COMPLETED**\nRuntime: ' + str(datetime.now() - start_time))
