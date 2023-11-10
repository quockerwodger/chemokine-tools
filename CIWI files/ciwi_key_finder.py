# CIWI created 2023 by Amy Scadden while a PhD student at the University of Otago

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

input_names = pd.read_csv('all_known.txt', sep = "\t", index_col = 'chemokine') #opens file used as argument
#print("Input table:\n", input_names)
#print(len(input_names))
input_sequences = input_names['sequence'].tolist()
#print(input_sequences)

# read in properties table

property_table_read = pd.read_csv('property_table.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')
#print("Table of properties:\n", property_table)

# read in consensus_table

class Consensus_type(Enum):
	POSITIVE = 'positive'
	GENERAL = 'general'
	NEGATIVE = 'negative'

while True:
	try:
		ask_table = Consensus_type(input("What consensus sequence would you like to use?:"))
	except ValueError:
		print("Invalid input. Please enter 'positive', 'general', or 'negative'.")
		continue
	else:
		break

if ask_table == Consensus_type.POSITIVE:
	consensus_table_read = pd.read_csv('positive_consensus_table.csv')
elif ask_table == Consensus_type.GENERAL:
	consensus_table_read = pd.read_csv('consensus_table.csv')
else:
	consensus_table_read = pd.read_csv('negative_consensus_table.csv')

consensus_table_read = consensus_table_read.rename(columns = {"Res_no":"Variable"})
consensus_table = consensus_table_read.set_index('Variable')
consensus_table.columns = [str(int(c)-1) for c in consensus_table.columns]
#print("Consensus table:\n", consensus_table)
consensus_sequence = consensus_table.loc['Res', :].values.flatten().tolist()
consensus_sequence = ''.join(consensus_sequence)
consensus_search = re.escape(consensus_sequence).replace('X','.')
#print("Consensus sequence:\n", consensus_sequence)
#print("Consensus search:\n", consensus_search)

# Create table for each input sequence and compare each to consensus table

all_res = pd.DataFrame() #creates empty dataframe
temp_res = pd.DataFrame()
comp_res = pd.DataFrame(index = consensus_table.index.values)
comp_scores_variable = pd.DataFrame(columns = input_names.index.values, index = consensus_table.index.values)
comp_scores_position = pd.DataFrame(columns = consensus_table.columns.values, index = input_names.index.values)
final_res = pd.DataFrame()
all_scores = pd.DataFrame()

chemokine = 0

for item in input_sequences:
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	#chemokine = 0
	if res == (len(item)):
		all_res = pd.concat([all_res, temp_res], axis = 0)
		temp_res = pd.DataFrame()
	#print(all_res)
	all_res.columns = consensus_table.columns
	for col in all_res:
		condition = np.logical_or(all_res[col] != consensus_table[col], consensus_table[col] == 'X')
		comp_res[col] = np.where(condition, 0, 1)
	#print(comp_res)
	col = chemokine
	row = chemokine
	means = comp_res.mean(axis = 1)
	#print(comp_res.mean(axis = 0))
	means_position = comp_res.mean(axis = 0)
	comp_scores_variable.iloc[:, col] = means
	comp_scores_position.iloc[col] = means_position
	#print(comp_scores_variable)
	#print(comp_scores_position)
	final_res = pd.concat([final_res, all_res], axis = 0)
	all_res = pd.DataFrame()
	all_scores = pd.concat([all_scores, comp_res], axis = 0)
	chemokine += 1

comp_scores_variable.loc['Mean'] = comp_scores_variable.mean()
comp_scores_position['Mean'] = comp_scores_position.mean(axis = 1)
#all_scores.columns = [str(int(c)+1) for c in all_scores.columns]
#all_res.columns = [str(int(c)+1) for c in all_res.columns]
#comp_scores_position.columns = [str(int(c)+1) for c in comp_scores_position.columns]
comp_scores_variable = comp_scores_variable.transpose()
print("All involvement:\n", all_scores)
print("Position scores:\n", comp_scores_position)
print("Variable scores:\n", comp_scores_variable)
comp_scores_variable.to_csv("Calculated means per variable.csv")
comp_scores_position.to_csv("Calculated means per consensus position.csv")
final_res.to_csv("All raw data.csv")
all_scores.to_csv("All raw scores.csv")
