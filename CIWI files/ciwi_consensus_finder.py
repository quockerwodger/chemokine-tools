# load required modules

import csv
import sys
import numpy as np
import scipy as sp
import pandas as pd
from collections import Counter
from collections import defaultdict

# read in input_sequences

input_files = open('binders.txt', 'r', encoding='UTF-8').readlines() #opens file used as argument
input_array = np.array(input_files) #converts file to array
input_sequences = [] #creates empty list
for sub in input_array: #takes each sub (value?) in the array and appends it to the empty list, removing the \n in the process
	input_sequences.append(sub.replace("\n", ""))
print("Input sequences:\n", input_sequences)
#print(len(input_sequences))

all_neg = pd.DataFrame()
temp_neg = pd.DataFrame()

# read in properties table

property_table_read = pd.read_csv('property_table.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')
print("Reference property table:\n", property_table)

#property_table = pd.read_csv('property_table.txt', sep = "\t", index_col = 'Variable') #reads property table as a dataframe
#print("Table of properties:\n", property_table)

# read in negative files
if len(sys.argv) > 2:
	negative_files = open(sys.argv[2], 'r', encoding = 'UTF-8').readlines()
	negative_array = np.array(negative_files)
	negative_sequences = []
	for sub in negative_array:
		negative_sequences.append(sub.replace("\n", ""))
	print("Negative sequences:\n", negative_sequences)
	for item in negative_sequences:
		res = 0
		while res < len (item):
			new_neg = {res:(property_table.loc[:, item[res]])}
			for key, value in new_neg.items():
				temp_neg[key] = value
			res = res + 1
		if res == (len(item)):
			all_neg = pd.concat([all_neg, temp_neg], axis = 0)
			temp_neg = pd.DataFrame()
	#print (all_neg)
else:
	print("Negative sequence:\nNo negative sequence provided")

# matching property table elements with sequences

all_res = pd.DataFrame() #creates empty dataframe
temp_res = pd.DataFrame()

for item in input_sequences: #reads through each input sequence, looks at each character (res) and finds the column in the property table with the same name. Places this information as a list into a holding dict (new_res) before extending the all_res dict with these values.
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	if res == (len(item)):
		all_res = pd.concat([all_res, temp_res], axis = 0)
		temp_res = pd.DataFrame()
print ("new_res:\n", new_res)
print ("temp_res:\n", temp_res)
print("all_res:\n", all_res)
print("len all_res:\n", len(all_res))
all_res.to_csv("full_table.csv")

# creating consensus properties (used ChatGPT for guidance, I highly recommend this service if you are struggling to get code to work)

# Count the number of occurrences of each value in each column
counts = all_res.apply(pd.Series.value_counts)

# Replace any value that appears in less than 60% of the sequences with an X
threshold = len(input_sequences) * 0.6
all_res = all_res.apply(lambda x: ['X' if counts[x.name][value] < threshold else value for value in x])

# Export frequency counts to new table
scores_df = all_res.apply(lambda x: [((counts[x.name][value])/len(input_sequences)) if value in counts[x.name] and counts[x.name][value] >= threshold else 0 for value in x])
print(threshold)
print(counts)

# Print the result
print("all_res:\n", all_res)
#print("scores_df:\n", scores_df)

all_res = all_res.sort_index(axis = 0)
scores_df = scores_df.sort_index(axis = 0)

for col in all_res:
    all_res = all_res.groupby('Variable').apply(lambda x: pd.concat([x[x[col] != "X"].sort_values(col), x[x[col] == "X"].sort_values(col)]))

for col in scores_df:
	scores_df = scores_df.groupby('Variable').apply(lambda x: x.apply(lambda col: col.sort_values(ascending = False)))

sorted_res = all_res.groupby('Variable').apply(lambda x: x.apply(lambda col: col.sort_values(ignore_index=False, key=lambda x: x!='X')))

print("sorted all_res:\n", sorted_res)
print("sorted scores_df:\n", scores_df)

# Creating the negative & positive consensus sequences
if all_neg.empty:
	print("Negative consensus:\nNo negative sequence provided")
else:
	neg_counts = all_neg.apply(pd.Series.value_counts)
	threshold_neg = len(negative_sequences)/2
	all_neg = all_neg.apply(lambda x: ['X' if neg_counts[x.name][value] < threshold_neg else value for value in x])
	#print(neg_counts)
	#print(threshold_neg)
	neg_scores = all_neg.apply(lambda x: [((neg_counts[x.name][value])/len(negative_sequences)) if value in neg_counts[x.name] and neg_counts[x.name][value] >= threshold_neg else 0 for value in x])
	all_neg = all_neg.sort_index(axis = 0)
	neg_scores = neg_scores.sort_index(axis = 0)
	#print("all_neg:\n", all_neg)
	#print("neg_scores:\n", neg_scores)
	for col in all_neg:
		all_neg = all_neg.groupby('Variable').apply(lambda x: pd.concat([x[x[col] != "X"].sort_values(col), x[x[col] == "X"].sort_values(col)]))
	for col in neg_scores:
		neg_scores = neg_scores.groupby('Variable').apply(lambda x: x.apply(lambda col: col.sort_values(ascending = False)))
	all_neg = all_neg.reset_index()
	neg_scores = neg_scores.reset_index()
	negative_consensus = all_neg.drop_duplicates(subset=['Variable'], keep = 'first')
	negative_final_scores = neg_scores.drop_duplicates(subset=['Variable'], keep = 'first')
	negative_consensus = negative_consensus.rename(columns = {"Variable":"Res_no"})
	negative_consensus = negative_consensus.set_index('Res_no')
	negative_final_scores = negative_final_scores.rename(columns = {"Variable":"Res_no"})
	negative_final_scores = negative_final_scores.set_index('Res_no')

sorted_res = sorted_res.reset_index()
scores_df = scores_df.reset_index()

final_consensus = sorted_res.drop_duplicates(subset=['Variable'], keep = 'last')
consensus_scores = scores_df.drop_duplicates(subset=['Variable'], keep = 'first')

print(final_consensus)

final_consensus = final_consensus.rename(columns = {"Variable":"Res_no"})
final_consensus = final_consensus.set_index('Res_no')
consensus_scores = consensus_scores.rename(columns = {"Variable":"Res_no"})
consensus_scores = consensus_scores.set_index('Res_no')
final_consensus.columns = [str(int(c)+1) for c in final_consensus.columns]
consensus_scores.columns = [str(int(c)+1) for c in consensus_scores.columns]

print("Consensus sequence:\n", final_consensus)
print("Consensus scores:\n", consensus_scores)

final_consensus.to_csv('consensus_table.csv')
consensus_scores.to_csv('consensus_scores.csv')

# Finding positive consensus
if all_neg.empty:
	print("Negative consensus:\nNo negative sequence provided")
else:
	#positive_consensus = pd.DataFrame(index = final_consensus.index.values)
	negative_consensus.columns = final_consensus.columns
	negative_final_scores.columns = final_consensus.columns
	#for col in final_consensus:
	print("Negative consensus:\n", negative_consensus)
	print("Negative scores:\n", negative_final_scores)
	negative_consensus.to_csv('negative_consensus_table.csv')
	negative_final_scores.to_csv('negative_consensus_scores.csv')
	positive_consensus = final_consensus.where(final_consensus != negative_consensus, 'X')
	print("Positive consensus:\n", positive_consensus)
	positive_scores = consensus_scores.where(final_consensus != negative_consensus, '0')
	print("Positive scores:\n", positive_scores)
	positive_consensus.to_csv('positive_consensus_table.csv')
	positive_scores.to_csv('positive_consensus_scores.csv')
