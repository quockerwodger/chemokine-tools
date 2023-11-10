# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 10:05:50 2023

@author: Josh
"""
import numpy as np
import pandas as pd
from collections import Counter
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
import pickle

#Loading in data

cluster_input = pd.read_csv("all_full_length.txt", sep = "\t", index_col = 'chemokine')
cluster_seq = cluster_input['sequence'].tolist()

full_input = pd.read_csv("all_chemokines.txt", sep = "\t", index_col = 'chemokine')
full_seq = full_input['sequence'].tolist()
    
property_table_read = pd.read_csv('property_table.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')

property_table2_read = pd.read_csv('property_table2.txt', sep = "\t") #reads property table as a dataframe
property_table2_read = property_table2_read.sort_values('Variable').reset_index(drop = True)
property_table2 = property_table2_read.set_index('Variable')

#Setting up data tables
#Models 1 & 2
model2_raw = pd.DataFrame()
temp_res = pd.DataFrame()
for item in cluster_seq:
	res = 0
	while res < len (item):
		new_res = {res:(property_table.loc[:, item[res]])}
		for key, value in new_res.items():
			temp_res[key] = value
		res = res + 1
	if res == (len(item)):
		model2_raw = pd.concat([model2_raw, temp_res], axis = 0)
		temp_res = pd.DataFrame()

hydropathy_consensus = [[-1.2833333333333334, 0.45613107278013365], [-1.5823529411764707, 1.8000768918601773], [-0.9199999999999999, 2.3576683396949623], [2.576470588235294, 0.3058823529411764], [2.5, 0.0], [3.3117647058823527, 0.6944166853352333], [-1.5299999999999998, 2.0935854412944317], [1.01, 2.506571363436517], [-0.725, 0.043301270189221974], [-1.7299999999999998, 1.7607100840286], [-4.147058823529412, 0.39722218861492076], [-3.435714285714286, 0.9961405113573822], [4.48, 0.07483314773547879], [-1.5529411764705883, 0.18823529411764706], [-0.6249999999999998, 3.4838018026288466], [-3.7714285714285714, 0.4043134770881402], [-4.027272727272727, 0.625359400797085], [4.235, 0.2815581645060218], [2.8, 0.0], [-0.33499999999999996, 3.1799803458512135], [-0.75625, 0.21785531322416726], [-4.0058823529411764, 0.22873202464968923], [-4.1625, 0.404467242184086], [-0.5090909090909091, 0.17814470856604933], [-4.175, 0.3381937314617171], [-3.566666666666667, 0.14907119849998596], [4.233333333333333, 0.2357022603955159], [2.5, 0.0], [1.8, 0.0]]
residue_consensus = ['X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'I', 'P', 'X', 'X', 'X', 'X', 'F', 'X', 'T', 'K', 'X', 'X', 'X', 'X', 'X', 'C', 'A']
volume_consensus = [[113.25, 1.4145081595145808], [127.00555555555555, 24.765353761123222], [112.56666666666665, 2.8365863678411443], [108.9470588235294, 1.7882352941176458], [108.5, 0.0], [190.9090909090909, 1.647838554235951], [92.25714285714285, 13.567382855621798], [192.0583333333333, 1.8241245997890427], [123.12999999999997, 28.372805642022783], [98.91999999999999, 12.2938575448609], [171.35333333333332, 2.546335056944049], [168.78181818181818, 2.727242424074074], [166.69999999999996, 0.01], [112.7, 0.0], [169.57272727272732, 3.7138783719965183], [133.53, 33.808802108326766], [172.43157894736842, 12.79444997405142], [166.7, 0.0], [189.9, 0.0], [143.66315789473683, 24.729318004649244], [115.88749999999999, 0.8230089610690741], [168.26470588235293, 0.7243180780573513], [170.94375, 2.853718790893742], [94.20999999999998, 37.098543098078665], [170.63529411764705, 2.669201505731766], [141.4857142857143, 2.6723069602491276], [155.27368421052628, 15.344460413255367], [108.5, 0.0], [86.92352941176469, 6.705882352941175]]
side_consensus = ['A', 'X', 'X', 'S', 'S', 'X', 'X', 'R', 'O', 'X', 'B', 'X', 'A', 'A', 'X', 'X', 'X', 'A', 'R', 'X', 'O', 'B', 'B', 'X', 'B', 'X', 'A', 'S', 'A']
polarity_consensus = [[0.0, 0.0], [-0.5833333333333334, 0.18633899812498245], ['X'], [0.0, 0.0], [0.0, 0], [0.0, 0.0], ['X'], ['X'], [-0.5, 0.0], ['X'], [0.9117647058823529, 0.19061002054140766], ['X'], [0.0, 0], [0.0, 0.0], ['X'], [0.6666666666666666, 0.23570226039551587], ['X'], [0.0, 0], [0.0, 0], ['X'], [-0.5, 0.0], [1.0, 0.0], [0.90625, 0.19515618744994995], ['X'], [0.9411764705882353, 0.16109486985446062], ['X'], [0.0, 0], [0.0, 0], [0.0, 0]]

#Setting up arrays
raw_group = model2_raw.groupby(model2_raw.index)
result = {}
for name, group in raw_group:
    result[name]=group

binder_hyd = result['Hydropathy']
binder_pol = result['Polarity']
binder_res = result['Residue']
binder_size = result['Size']
binder_side = result['Side_chain']

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

#Scoring side chain type
binder_side_scores = []
binary_match(binder_side, side_consensus, binder_side_scores)

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

#Scoring size
binder_size_scores = []
continuous_match(binder_size, volume_consensus, binder_size_scores)

#Scoring polarity
binder_pol_scores = []
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

#Side chains
binder_side_averages = []
for value in binder_side_scores:
    average = np.mean(value)
    binder_side_averages.extend([average])
    
#Hydropathy
def cont_mean(test, output, threshold, weight):
    for data in test:
        weights = [1 if x <= threshold else weight for x in data]
        weighted_average = np.average(data, weights=weights)
        output.extend([weighted_average])
binder_hyd_averages = []
cont_mean(binder_hyd_scores, binder_hyd_averages, 5, 0.8)

#Size
binder_size_averages = []
non_size_averages = []
cont_mean(binder_size_scores, binder_size_averages, 100, 0.05)

#Polarity
binder_pol_averages = []
non_pol_averages = []
cont_mean(binder_pol_scores, binder_pol_averages, 25, 0.2)

model1_table = pd.concat([pd.DataFrame(binder_hyd_averages, index = cluster_input.index.values, columns = ['Hyd']), 
                          pd.DataFrame(binder_side_averages, index = cluster_input.index.values, columns = ['Side']), 
                          pd.DataFrame(binder_size_averages, index = cluster_input.index.values, columns = ['Size']), 
                          pd.DataFrame(binder_res_averages, index = cluster_input.index.values, columns = ['Res']), 
                          pd.DataFrame(binder_pol_averages, index = cluster_input.index.values, columns = ['Pol'])], axis = 1, ignore_index = True)
model1_table = model1_table.reset_index()
model2_table = pd.concat([pd.DataFrame(data = binder_hyd_scores), pd.DataFrame(data = binder_side_scores), pd.DataFrame(binder_res_scores), pd.DataFrame(binder_size_scores), pd.DataFrame(binder_pol_scores)], axis = 1)

model1_features = []
for index, row in model1_table.iterrows():
    row_array = row.to_numpy()
    new_array = (row_array[1:]).astype(float)
    model1_features.append(new_array)
model2_features = []
for index, row in model2_table.iterrows():
    row_array = row.to_numpy()
    new_array = (row_array).astype(float)
    model2_features.append(new_array)

### Model 3 ###
        
model3_array = []
temp_res = pd.DataFrame()
for item in cluster_seq:
    res = 0
    temp_array = []
    while res < len (item):
        new_res = {res:(property_table2.loc[:, item[res]])}
        for key, value in new_res.items():
            temp_res[key] = value
            res = res + 1
    if res == (len(item)):
        for index, row in temp_res.iterrows():
            row_array = (row.to_numpy()).astype(float)
            temp_array.extend(row_array)
        model3_array.append(temp_array)
        temp_res = pd.DataFrame()

### Model 4 ###

model4_array = []
for item in full_seq:
    res = 0
    temp_array = []
    while res < len (item):
        new_res = {res:(property_table2.loc[:, item[res]])}
        for key, value in new_res.items():
            temp_res[key] = value
            res = res + 1
    if res == (len(item)):
        for index, row in temp_res.iterrows():
            row_array = (row.to_numpy()).astype(float)
            temp_array.extend(row_array)
        model4_array.append(temp_array)
        temp_res = pd.DataFrame()

### Running models ###

model1 = pickle.load(open('model1.sav', 'rb'))
model2 = pickle.load(open('model2.sav', 'rb'))
model3 = pickle.load(open('model3.sav', 'rb'))
model4 = pickle.load(open('model4.sav', 'rb'))

model1_out = model1.predict(model1_features)
model1_prob = model1.predict_proba(model1_features)
print("Model 1")
for out, prob in zip(model1_out, model1_prob):
    print(out, ",", prob)
    
model2_out = model2.predict(model2_features)
model2_prob = model2.predict_proba(model2_features)
print("Model 2")
for out, prob in zip(model2_out, model2_prob):
    print(out, ",", prob)

model3_out = model3.predict(model3_array)
model3_prob = model3.predict_proba(model3_array)
print("Model 3")
for out, prob in zip(model3_out, model3_prob):
    print(out, ",", prob)
    
model4_out = model4.predict(model4_array)
model4_prob = model4.predict_proba(model4_array)
print("Model 4")
for out, prob in zip(model4_out, model4_prob):
    print(out, ",", prob)