# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 20:12:38 2023

@author: Josh
"""
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error

input_files = open("chemokines.txt").readlines() #opens file used as argument
input_array = np.array(input_files) #converts file to array
input_sequences = [] #creates empty list
for sub in input_array: #takes each sub (value?) in the array and appends it to the empty list, removing the \n in the process
	input_sequences.append(sub.replace("\n", ""))
#non_input = pd.read_csv("non_binders.txt", sep = "\t", index_col = 'chemokine')
#non_seq = non_input['sequence'].tolist()

property_table_read = pd.read_csv('property_table2.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')

chemokine_array = []
temp_res = pd.DataFrame()
for item in input_sequences:
    res = 0
    temp_array = []
    while res < len (item):
        new_res = {res:(property_table.loc[:, item[res]])}
        for key, value in new_res.items():
            temp_res[key] = value
            res = res + 1
    if res == (len(item)):
        for index, row in temp_res.iterrows():
            row_array = (row.to_numpy()).astype(float)
            temp_array.extend(row_array)
        chemokine_array.append(temp_array)
        temp_res = pd.DataFrame()
binding_status = []
for i, item in enumerate(chemokine_array):
    if i< 20:
        value = 'Binder'
    else:
        value = 'Non-binder'
    binding_status.extend([value])

kd_all_files = open("all_kd.txt").readlines() #opens file used as argument
kd_all_array = np.array(kd_all_files) #converts file to array
all_kd = [] #creates empty list
for sub in kd_all_array: #takes each sub (value?) in the array and appends it to the empty list, removing the \n in the process
	all_kd.append(float(sub.replace("\n", "")))
#Machine learning
train4, test4, train_labels4, test_labels4 = train_test_split(chemokine_array, binding_status, test_size = 0.25, random_state = 0)

clf4 = MLPClassifier(max_iter = int(1e19), activation = 'logistic', solver = 'adam', alpha = 0.01, learning_rate_init = 0.0001)
print(clf4)

model = clf4.fit(train4, train_labels4)

preds = clf4.predict(test4)
print(preds)

print(accuracy_score(test_labels4, preds))

#scores = cross_val_score(clf4, chemokine_array, binding_status, cv=20)
#print(scores.mean(), scores.std())

import pickle
pickle.dump(model, open('model4.sav', 'wb'))
