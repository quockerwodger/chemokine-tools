# -*- coding: utf-8 -*-

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

cluster_input = pd.read_csv("all_chemokines.txt", sep = "\t", index_col = 'chemokine')
cluster_seq = cluster_input['sequence'].tolist()

full_input = pd.read_csv("all_long2.txt", sep = "\t", index_col = 'chemokine')
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

hydropathy_consensus = [[-1.2833333333333334, 0.45613107278013365], [-0.9000000000000001, 2.7011461412131106], [-0.7571428571428571, 2.413369789629554], [2.6666666666666665, 0.47609522856952335], [2.5, 0.0], [3.3388888888888886, 0.6840556254906703], [-1.6428571428571426, 2.1045464229278723], [0.9285714285714286, 2.4731208079291727], [-0.725, 0.043301270189221974], [-1.7238095238095237, 1.7185001106736906], [-4.166666666666667, 0.39440531887330776], [-3.435714285714286, 0.9961405113573822], [4.4625, 0.09921567416492208], [-1.5529411764705883, 0.18823529411764706], [-0.4619047619047617, 3.477201418557087], [-3.78, 0.39191835884530846], [-3.983333333333334, 0.6162160515778716], [4.247619047619048, 0.2805081233337742], [4.5, 0.0], [2.8, 0.0], [-0.35238095238095235, 3.1043163571632766], [-0.75625, 0.21785531322416726], [-4.0, 0.223606797749979], [-4.1625, 0.404467242184086], [-0.5000000000000001, 0.17320508075688773], [-4.135294117647058, 0.364516079643051], [-3.5842105263157897, 0.163072982998207], [4.231578947368422, 0.22953644723533187], [2.5, 0.0], [1.8, 0.0]]
residue_consensus = ['X', 'X', 'X', 'C', 'C', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'I', 'P', 'X', 'X', 'X', 'X', 'I', 'F', 'X', 'T', 'K', 'X', 'X', 'X', 'X', 'X', 'C', 'A']
volume_consensus = [[113.25, 1.4145081595145808], [134.41904761904763, 29.39570132240134], [112.3125, 2.917592115083942], [108.9470588235294, 1.7882352941176458], [108.5, 0.0], [190.9090909090909, 1.647838554235951], [95.69333333333333, 18.360445407330282], [192.0583333333333, 1.8241245997890427], [125.20476190476187, 29.20230708743788], [106.0190476190476, 22.36808140542316], [171.48125, 2.514761109429687], [168.78181818181818, 2.727242424074074], [166.69999999999996, 0.01], [112.60588235294118, 0.37647058823529617], [169.57272727272732, 3.7138783719965183], [135.2, 33.828728281489646], [172.43157894736842, 12.79444997405142], [166.69999999999996, 0.01], [166.69999999999996, 0.01], [189.9, 0.0], [142.28499999999997, 24.840476545348317], [115.88749999999999, 0.8230089610690741], [168.28333333333333, 0.7080881928749354], [170.94375, 2.853718790893742], [92.58571428571426, 36.92600433188469], [170.63529411764705, 2.669201505731766], [141.4857142857143, 2.6723069602491276], [151.37142857142854, 20.50839881045586], [108.5, 0.0], [86.92352941176469, 6.705882352941175]]
side_consensus = ['A', 'X', 'X', 'S', 'S', 'X', 'X', 'X', 'X', 'X', 'B', 'X', 'A', 'A', 'X', 'X', 'X', 'A', 'A', 'R', 'X', 'O', 'B', 'B', 'X', 'B', 'X', 'A', 'S', 'A']
polarity_consensus = [[0.0, 0.0], ['X'], ['X'], [0.0, 0.0], [0.0, 0], [0.0, 0.0], ['X'], ['X'], ['X'], ['X'], [0.9166666666666666, 0.18633899812498245], ['X'], [0.0, 0], [0.0, 0.0], ['X'], [0.7, 0.24494897427831783], ['X'], [0.0, 0], [0.0, 0], [0.0, 0], ['X'], [-0.5, 0.0], [1.0, 0.0], [0.90625, 0.19515618744994995], ['X'], [0.9411764705882353, 0.16109486985446062], ['X'], [0.0, 0], [0.0, 0], [0.0, 0]]

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