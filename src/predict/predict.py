import pandas
from sklearn.cluster import KMeans
from yellowbrick.cluster import KElbowVisualizer
import sys
import os 
import numpy as np

import warnings
warnings.filterwarnings('ignore')

query_output = sys.argv[1]
output_file = sys.argv[2]

def get_idx_min(df, column_name, cov_column_name):
	df1 = df.groupby([column_name])[cov_column_name].agg(['mean','count'])
	mean_cov = []

	for cluster_i in range(len(df1)):      
		mean_cov.append(df1.loc[(cluster_i,)]['mean']) 
		
	idx_min = 0
	for i in range(1, len(mean_cov)):
		if mean_cov[i] < mean_cov[idx_min]:
			idx_min = i

	return idx_min

def predict(df, idx_min, nclusters, output_file):	
	df['predict'] = np.where(df['KM{}'.format(nclusters)] == idx_min, 0, 1)
	with open(output_file, 'w') as f:
		for i in range(len(df['predict'])):
			if df['predict'][i] == 1:
				f.write(df['name'][i]+'\n')
				# print(df['species'][i], df['predict'][i])

def get_all_names(df, output_file):
	with open(output_file, 'w') as f:
		for i in range(len(df['name'])):
			f.write(df['name'][i]+'\n')

km = KMeans()
df = pandas.read_csv(query_output)
X = df[['coverage_of_uniq_sigs']]
km_visualizer = KElbowVisualizer(km, k=(2,12), metric='distortion', timings=False)
try:
	km_a = km_visualizer.fit(X)
	km_optimal_n_clusters = km_a.elbow_value_
	print("Optimal number of clusters = ", km_optimal_n_clusters)

	km = KMeans(n_clusters=km_optimal_n_clusters)
	df['KM{}'.format(km_optimal_n_clusters)] = km.fit_predict(X)
	idx_min = get_idx_min(df, 'KM{}'.format(km_optimal_n_clusters), 'coverage_of_uniq_sigs')
	print("Lowest coverage group id:", idx_min)
	predict(df, idx_min, km_optimal_n_clusters, output_file)
except ValueError:
	get_all_names(df, output_file)



