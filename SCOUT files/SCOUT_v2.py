import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score

binder_input = pd.read_csv("binders.txt", sep = "\t", index_col = 'chemokine')
binder_seq = binder_input['sequence'].tolist()
non_input = pd.read_csv("non_binders.txt", sep = "\t", index_col = 'chemokine')
non_seq = non_input['sequence'].tolist()

property_table_read = pd.read_csv('property_table2.txt', sep = "\t") #reads property table as a dataframe
property_table_read = property_table_read.sort_values('Variable').reset_index(drop = True)
property_table = property_table_read.set_index('Variable')

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
        
binder_group = binder_table.groupby(binder_table.index)
result = {}
for name, group in binder_group:
    result[name]=group

binder_hyd = result['Hydropathy'].reset_index().drop(columns = ['Variable'])
binder_pol = result['Polarity'].reset_index().drop(columns = ['Variable'])
binder_size = result['Size'].reset_index().drop(columns = ['Variable'])
all_binder = pd.concat([binder_hyd, binder_pol, binder_size], axis = 1)
non_group = non_table.groupby(non_table.index)
result = {}
for name, group in non_group:
    result[name]=group

non_hyd = result['Hydropathy'].reset_index().drop(columns = ['Variable'])
non_pol = result['Polarity'].reset_index().drop(columns = ['Variable'])
non_size = result['Size'].reset_index().drop(columns = ['Variable'])
all_non = pd.concat([non_hyd, non_pol, non_size], axis = 1)

all_data = pd.concat([all_binder, all_non], axis = 0)
# Preparing data
all_data = all_data.reset_index()
binding_types = ['Binder', 'Non-binder']
binding_status = []
for index, row in all_data.iterrows():
    if index < 20:
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
    print(row_array)
    new_array = (row_array[1:]).astype(float)
    features.append(new_array)
print(features)

kd_all_files = open("all_kd.txt").readlines() #opens file used as argument
kd_all_array = np.array(kd_all_files) #converts file to array
all_kd = [] #creates empty list
for sub in kd_all_array: #takes each sub (value?) in the array and appends it to the empty list, removing the \n in the process
	all_kd.append(float(sub.replace("\n", "")))
print(all_kd)
# Training & testing the model 
train3, test3, train_labels3, test_labels3 = train_test_split(features, binding_status, test_size = 0.25, random_state = 0)

clf3 = MLPClassifier(hidden_layer_sizes = (1000,), max_iter = int(1e19), activation = 'logistic', solver = 'adam', alpha = 0.01, learning_rate_init = 0.0001)
print(clf3)

model = clf3.fit(train3, train_labels3)

preds = clf3.predict(test3)
print(preds)
print(accuracy_score(test_labels3, preds))

#scores = cross_val_score(clf3, features, binding_status, cv=20)
#print(scores.mean(), scores.std())
import pickle
pickle.dump(model, open('model3.sav', 'wb'))
