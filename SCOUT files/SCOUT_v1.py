# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from collections import Counter
import sklearn

###Input files###
binder_input = pd.read_csv("binders.txt", sep = "\t", index_col = 'chemokine')
binder_seq = binder_input['sequence'].tolist()
non_input = pd.read_csv("non_binders.txt", sep = "\t", index_col = 'chemokine')
non_seq = non_input['sequence'].tolist()

property_table_read = pd.read_csv('property_table.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')

###Converting sequences into a variable table###
binder_table = pd.DataFrame()
temp_res = pd.DataFrame()
for item in binder_seq:
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	if res == (len(item)):
		binder_table = pd.concat([binder_table, temp_res], axis = 0)
		temp_res = pd.DataFrame()

non_table = pd.DataFrame()
temp_res = pd.DataFrame()
for item in non_seq:
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	if res == (len(item)):
		non_table = pd.concat([non_table, temp_res], axis = 0)
		temp_res = pd.DataFrame()

#Splitting variables into separate tables
grouped = binder_table.groupby(binder_table.index)
result = {}
for name, group in grouped:
    result[name]=group

hydropathy_df = result['Hydropathy']
polarity_df = result['Polarity']
residue_df = result['Residue']
size_df = result['Size']
side_df = result['Side_chain']

###Scoring each table###
#Residue
residue_consensus = []
for col in residue_df:
    column = residue_df[col].to_numpy()
    occurence_count = (Counter(column))
    top_count = occurence_count.most_common(1)
    if top_count[0][1]>=(len(column)*0.6):
        residue_consensus.extend(top_count[0][0])
    else:
        residue_consensus.extend("X")

#Volume
def find_mean(data, cutoff):
    med_val = np.median(data)
    diff_from_median = []
    for item in data:
        abs_diff = [abs(item-med_val)]
        diff_from_median.extend(abs_diff)
    med_dev = np.median(diff_from_median)
    mean_dev = np.mean(diff_from_median)
    values_using = []
    for item in data:
        if med_dev ==0:
            if mean_dev == 0:
                z_score = 0
            else:
                z_score = (item-med_val)/(1.253314*mean_dev)
        else:
            z_score = (item-med_val)/(1.486*med_dev)
        if -(cutoff) <= z_score <= cutoff:
            values_using.extend([item])
    mean_value = np.mean(values_using)
    if 0.0<np.std(values_using)<0.0001:
        std_dev = 0.01
    else:
        std_dev = np.std(values_using)
    output = [mean_value, std_dev]
    return(output)
volume_consensus = []
for col in size_df:
    column = size_df[col].to_numpy()
    column_float = column.astype(float)
    col_mean = find_mean(column_float, 2)
    volume_consensus.append(col_mean)

#Hydropathy
hydropathy_consensus= []
for col in hydropathy_df:
    column = (hydropathy_df[col].to_numpy()).astype(float)
    col_mean = find_mean(column, 2)
    hydropathy_consensus.append(col_mean)

#Polarity
polarity_consensus = []
for col in polarity_df:
    column = (polarity_df[col].to_numpy()).astype(float)
    occurence_count = (Counter(column))
    top_count = occurence_count.most_common(1)
    #std_dev = np.std(column)
    if len(occurence_count) == 1:
        final_value = top_count[0][0]
        polarity_consensus.append([final_value, 0])
    if len(occurence_count)>1:
        counts = occurence_count.most_common(len(occurence_count))
        positive_polar = [item for item in counts if 0.5 in item]
        positive_charged = [item for item in counts if 1.0 in item]
        negative_polar = [item for item in counts if -0.5 in item]
        negative_charged = [item for item in counts if -1.0 in item]
        non_polar = [item for item in counts if 0.0 in item]
        if len(positive_polar) > 0 and len(positive_charged) > 0 and (positive_polar[0][1]+ positive_charged[0][1]) >= len(column)*0.6:
            final_value = ((positive_polar[0][0]*positive_polar[0][1])+(positive_charged[0][0]*positive_charged[0][1]))/(positive_polar[0][1]+positive_charged[0][1])
            new_set = [positive_polar[0], positive_charged[0]]
            new_array = [value for value, count in new_set for _ in range(count)]
            std_dev = np.std(new_array)
            polarity_consensus.append([final_value, std_dev])
        elif len(negative_polar) > 0 and len(negative_charged) > 0 and (negative_polar[0][1]+ negative_charged[0][1]) >= len(column)*0.6:
            final_value = ((negative_polar[0][0]*negative_polar[0][1])+(negative_charged[0][0]*negative_charged[0][1]))/(negative_polar[0][1]+negative_charged[0][1])
            new_set = [negative_polar[0], negative_charged[0]]
            new_array = [value for value, count in new_set for _ in range(count)]
            std_dev = np.std(new_array)
            polarity_consensus.append([final_value, std_dev])
        else:
            if top_count[0][1]>= len(column)*0.6:
                final_value = top_count[0][0]
                new_array = [value for value, count in top_count for _ in range(count)]
                std_dev = np.std(new_array)
                polarity_consensus.append([final_value, std_dev])
            else:
                final_value = "X"
                polarity_consensus.append([final_value])
#Side chain
side_consensus = []
for col in side_df:
    column = side_df[col].to_numpy()
    occurence_count = (Counter(column))
    top_count = occurence_count.most_common(1)
    if top_count[0][1]>=(len(column)*0.6):
       side_consensus.extend(top_count[0][0])
    else:
        side_consensus.extend("X")


### Scoring inputs ###
#Setting up arrays
binder_group = binder_table.groupby(binder_table.index)
result = {}
for name, group in binder_group:
    result[name]=group

binder_hyd = result['Hydropathy']
binder_pol = result['Polarity']
binder_res = result['Residue']
binder_size = result['Size']
binder_side = result['Side_chain']

non_group = non_table.groupby(non_table.index)
result = {}
for name, group in non_group:
    result[name]=group

non_hyd = result['Hydropathy']
non_pol = result['Polarity']
non_res = result['Residue']
non_size = result['Size']
non_side = result['Side_chain']

#Scoring residue id
def binary_match(test, consensus, output):
    for index, row in test.iterrows():
        row_array = row.to_numpy()
        score_array = []
        for item in zip(row_array, consensus):
            if item[0] == item [1]:
                score_array.extend([0])
            elif item [1] == 'X':
                score_array.extend([0])
            else:
                score_array.extend([1])
        output.append(score_array)
    return(output)
binder_res_scores = []
binary_match(binder_res, residue_consensus, binder_res_scores)
non_res_scores = []
binary_match(non_res, residue_consensus, non_res_scores)

#Scoring side chain type
binder_side_scores = []
non_side_scores = []
binary_match(binder_side, side_consensus, binder_side_scores)
binary_match(non_side, side_consensus, non_side_scores)

#Scoring hydropathy
def continuous_match(test, consensus, output):
    for index, row in test.iterrows():
        row_array = row.to_numpy().astype(float)
        score_array = []
        for item in zip(row_array, consensus):
            if item[1][1] == 0:
                if item[0] == item[1][0]:
                    score_array.extend([0])
                else:
                    if (item[0]-item[1][0])>0:
                        deviation = 0.01
                    else:
                        deviation = 0.01
                    score = abs(item[0]-item[1][0])/deviation
                    score_array.extend([score])
            else:
                score = (abs(item[0]-item[1][0])/item[1][1])
                score_array.extend([score])
            if 0<item[1][1]<0.01:
                print(item)
                print(score)
        output.append(score_array)
binder_hyd_scores = []
continuous_match(binder_hyd, hydropathy_consensus, binder_hyd_scores)
non_hyd_scores = []
continuous_match(non_hyd, hydropathy_consensus, non_hyd_scores)

#Scoring size
binder_size_scores = []
non_size_scores = []
continuous_match(binder_size, volume_consensus, binder_size_scores)
continuous_match(non_size, volume_consensus, non_size_scores)

#Scoring polarity
binder_pol_scores = []
non_pol_scores = []
def polar_match(test, consensus, output):
    for index, row in test.iterrows():
        row_array = row.to_numpy().astype(float)
        score_array = []
        for item in zip(row_array, consensus):
            if len(item[1])>1:
                if item[1][1] == 0:
                    if item[0] == item[1][0]:
                        score_array.extend([0])
                    else:
                        if (item[0]-item[1][0])>0:
                            deviation = 0.01
                        else:
                            deviation = 0.01
                        score = abs(item[0]-item[1][0])/deviation
                        score_array.extend([score])
                else:
                    score = (abs(item[0]-item[1][0])/item[1][1])
                    score_array.extend([score])
            if len(item[1])<2:
                if item[1][0] == 'X':
                    score_array.extend([0])
        output.append(score_array)
polar_match(binder_pol, polarity_consensus, binder_pol_scores)
polar_match(non_pol, polarity_consensus, non_pol_scores)

#Combining outputs for future reference
'''
all_hyd = pd.concat([pd.DataFrame(binder_hyd_scores), pd.DataFrame(non_hyd_scores)])
all_pol = pd.concat([pd.DataFrame(binder_pol_scores), pd.DataFrame(non_pol_scores)])
all_size = pd.concat([pd.DataFrame(binder_size_scores), pd.DataFrame(non_size_scores)])
all_res = pd.concat([pd.DataFrame(binder_res_scores), pd.DataFrame(non_res_scores)])
all_side = pd.concat([pd.DataFrame(binder_side_scores), pd.DataFrame(non_side_scores)]) 
'''
###Creating average scores###
#Residues
binder_res_averages = []
for value in binder_res_scores:
    average = np.mean(value)
    binder_res_averages.extend([average])
non_res_averages = []
for value in non_res_scores:
    average = np.mean(value)
    non_res_averages.extend([average])

#Side chains
binder_side_averages = []
for value in binder_side_scores:
    average = np.mean(value)
    binder_side_averages.extend([average])
non_side_averages = []
for value in non_side_scores:
    average = np.mean(value)
    non_side_averages.extend([average])
    
#Hydropathy
def cont_mean(test, output, threshold, weight):
    for data in test:
        weights = [1 if x <= threshold else weight for x in data]
        weighted_average = np.average(data, weights=weights)
        output.extend([weighted_average])
binder_hyd_averages = []
non_hyd_averages = []
cont_mean(binder_hyd_scores, binder_hyd_averages, 5, 0.8)
cont_mean(non_hyd_scores, non_hyd_averages, 5, 0.8)

#Size
binder_size_averages = []
non_size_averages = []
cont_mean(binder_size_scores, binder_size_averages, 100, 0.05)
cont_mean(non_size_scores, non_size_averages, 100, 0.05)

#Polarity
binder_pol_averages = []
non_pol_averages = []
cont_mean(binder_pol_scores, binder_pol_averages, 25, 0.2)
cont_mean(non_pol_scores, non_pol_averages, 25, 0.2)

###Creating a table with all averages & chemokine names, plus binding status###
all_hyd = pd.concat([pd.DataFrame(data = binder_hyd_averages, index=binder_input.index.values, columns=['Hyd']), pd.DataFrame(data = non_hyd_averages, index=non_input.index.values, columns=['Hyd'])])
all_side = pd.concat([pd.DataFrame(data = binder_side_averages, index=binder_input.index.values, columns=['Side']), pd.DataFrame(data = non_side_averages, index=non_input.index.values, columns=['Side'])])
all_res = pd.concat([pd.DataFrame(data = binder_res_averages, index=binder_input.index.values, columns=['Res']), pd.DataFrame(data = non_res_averages, index=non_input.index.values, columns=['Res'])])
all_size = pd.concat([pd.DataFrame(data = binder_size_averages, index=binder_input.index.values, columns=['Size']), pd.DataFrame(data = non_size_averages, index=non_input.index.values, columns=['Size'])])
all_pol = pd.concat([pd.DataFrame(data = binder_pol_averages, index=binder_input.index.values, columns=['Pol']), pd.DataFrame(data = non_pol_averages, index=non_input.index.values, columns=['Pol'])])
all_averages = pd.concat([all_hyd, all_side, all_res, all_size, all_pol], axis = 1, ignore_index = True)
all_averages = all_averages.reset_index()

### Creating a table with all raw_data ###
all_binder = pd.concat([pd.DataFrame(data = binder_hyd_scores), pd.DataFrame(data = binder_side_scores), pd.DataFrame(binder_res_scores), pd.DataFrame(binder_size_scores), pd.DataFrame(binder_pol_scores)], axis = 1)
all_non = pd.concat([pd.DataFrame(data = non_hyd_scores), pd.DataFrame(data = non_side_scores), pd.DataFrame(non_res_scores), pd.DataFrame(non_size_scores), pd.DataFrame(non_pol_scores)], axis = 1)
all_data = pd.concat([all_binder, all_non], ignore_index = True)

### Machine learning model ###
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import mean_squared_error

# Preparing data
binding_types = ['Binder', 'Non-binder']
binding_status = []
for index, row in all_data.iterrows():
    if index < 21:
        value = 'Binder'
    else:
        value = 'Non-binder'
    binding_status.extend([value])
feature_names = []
for col in all_data:
    if col != 'index':
        feature_names.extend([col])
features = []
for index, row in all_data.iterrows():
    row_array = row.to_numpy()
    new_array = (row_array).astype(float)
    features.append(new_array)

mean_features = []
for index, row in all_averages.iterrows():
    row_array = row.to_numpy()
    new_array = (row_array[1:]).astype(float)
    mean_features.append(new_array)
print(mean_features)

#Model 1
clf1 = MLPClassifier(hidden_layer_sizes = (1000,), max_iter = int(1e19), activation = 'logistic', solver = 'adam', alpha = 0.01, learning_rate_init = 0.0001)
train1, test1, train_labels1, test_labels1 = train_test_split(mean_features, binding_status, test_size = 0.25, random_state = 42, stratify = binding_status)
model1 = clf1.fit(train1, train_labels1)
preds1 = clf1.predict(test1)
print(preds1)

print(accuracy_score(test_labels1, preds1))

scores1 = cross_val_score(clf1, mean_features, binding_status, cv = StratifiedKFold(10, shuffle = True, random_state = 42))
print(scores1.mean(), scores1.std())

probs1 = model1.predict_proba(test1)

# Model 2
train2, test2, train_labels2, test_labels2 = train_test_split(features, binding_status, test_size = 0.25, random_state = 0)

clf2 = MLPClassifier(hidden_layer_sizes = (1000,), max_iter = int(1e19), activation = 'logistic', solver = 'adam', alpha = 0.01, learning_rate_init = 0.0001)
model2 = clf2.fit(train2, train_labels2)

preds2 = clf2.predict(test2)
print(preds2)

print(accuracy_score(test_labels2, preds2))

scores2 = cross_val_score(clf2, features, binding_status, cv = StratifiedKFold(10, shuffle = True, random_state = 42))
print(scores2.mean(), scores2.std())

#Saving models
import pickle
filename1 = 'model1.sav'
filename2 = 'model2.sav'

pickle.dump(model1, open(filename1, 'wb'))
pickle.dump(model2, open(filename2, 'wb'))


