import numpy as np
import pandas as pd
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from sklearn.cluster import KMeans
from hdbscan import HDBSCAN

def calc_kernel_PCA(matrix, kernel = "precomputed") :
	from sklearn.decomposition import KernelPCA
	pca = KernelPCA(kernel = kernel)
	gram = 1 - (matrix / matrix.max())

	xy = pca.fit_transform(gram)

	return(xy)

def calculate_linkage(mtx, index_vals) :
	mtx_df = pd.DataFrame(mtx, index = index_vals, columns = index_vals)
	mtx_corr = mtx_df.corr()
	mtx_dism = 1 - mtx_corr
	linkage = hc.linkage(sp.distance.squareform(mtx_dism), method='average')

	return(linkage)

def cut_dist_tree(linkage, n_clusts) :
	hc.cut_tree(linkage, n_clusters = n_clust)

def kmeans(data, k = 30) :
	kmeans = KMeans(n_clusters = k, random_state=0).fit(data)

	return(kmeans.labels_)

def hdbscan(data, **kwargs) :
	clusterer = HDBSCAN(**kwargs)
	clusterer.fit(data)
	return(clusterer.labels_)

def tcrdiv(data, by_factor = None) :
	'''Calculate the TCRDiv as desribed by Dash et al.  
	Derived from TCRDist github: analyze_overlap_compute_simpsons.py'''
	div_dict = {}
	if by_factor :
		if by_factor not in data.obs.columns :
			raise ValueError("{by_fac} not in .obs".format(by_fac = by_factor))
		unique_factors = set(data.obs[by_factor])

		for fac in unique_factors :
			data_sub = data[data.obs[by_factor] == fac]
			div = _calc_tcrdiv(data_sub.X)
			div_dict[fac] = div

	else :
		div = _calc_tcrdiv(data.X)
		div_dict["all"] = div

	return(div_dict)

# Modified such that it dynamically calculates the standard deviationof the input matrix
def _calc_tcrdiv(matrix) :
	overlap_sum = 0
	total_sum = 0
	
	sd = np.std(matrix)
	for i in range(matrix.shape[0]) :
		for j in range(i + 1, matrix.shape[0]) :
			dist = matrix[i,j]
			if dist == 0 :
				overlap_sum += 1
			else :
				overlap_sum += np.exp(-1 * (dist/sd)**2)
			total_sum += 1

	p0 = overlap_sum / total_sum
	diversity = 1 / p0
	return(diversity.item())

