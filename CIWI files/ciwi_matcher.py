# load required modules

import csv
import sys
import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr
import re
import Bio
from Bio import motifs
import itertools
from enum import Enum

# read in input_sequences

input_names = pd.read_csv('all_chemokines.txt', sep = "\t", index_col = 'chemokine', encoding = "utf-8") #opens file used as argument
print("Input table:\n", input_names)
#print(len(input_names))
input_sequences = input_names['sequence'].tolist()
print("Input sequences:\n", input_sequences)

# read in properties table

property_table_read = pd.read_csv('property_table.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')
print("Table of properties:\n", property_table)

# read in consensus_table

input_consensus = pd.read_csv('consensus_table.csv', index_col = 'Variable')
print("Input consensus:\n", input_consensus)

# calculating number of values in input_consensus
match_count = 0
for col in input_consensus:
	match_count += ((input_consensus[col] != 'X').sum())
print("Matches in consensus:\n", match_count)

# Calculate match scores
all_res = pd.DataFrame()
temp_res = pd.DataFrame()
calc_table = pd.DataFrame()
item_count = 0
item_score = 0
chemokine = 0
match_list = []
all_raw = pd.DataFrame()
all_scores = pd.DataFrame()

for item in input_sequences:
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	if res == (len(item)):
		all_res = pd.concat([all_res, temp_res], axis = 0)
		temp_res = pd.DataFrame()
	all_res.columns = input_consensus.columns
	for col in all_res:
		calc_table[col] = np.where(all_res[col] == input_consensus[col], 1, 0)
	#print(all_res)
	#print(calc_table)
	for col in calc_table:
		item_count += ((calc_table[col] == 1).sum())
	item_score = item_count/match_count
	match_list.append(item_score)
	item_count = 0
	item_score = 0
	all_raw = pd.concat([all_raw, all_res], axis = 0)
	all_scores = pd.concat([all_scores, calc_table], axis = 0)
	all_res = pd.DataFrame()
	calc_table = pd.DataFrame()
	chemokine += 1

print("Match list:\n", match_list)

match_scores = pd.DataFrame(index = input_names.index.values, columns = ['Match Score'])
match_scores['Match Score'] = match_list
print("Match scores:\n", match_scores)
match_scores.to_csv("Match scores.csv")
all_raw.to_csv("All properties.csv")
all_scores.to_csv("All matches.csv")
