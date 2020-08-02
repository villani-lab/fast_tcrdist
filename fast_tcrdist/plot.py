import seaborn as sns
import pandas as pd
from . import useful
from .tcrdist import stats as stat
import matplotlib.pyplot as plt
import numpy as np	
import logomaker
from collections import Counter

def plot_distance_matrix(
	data,
	cell_id,
	row_colors : str = None,
	**kwargs) :
	
	if row_colors :
		if row_colors not in  data.obs.columns :
			raise ValueError("{row_col} not in .obs".format(row_col = row_colors))
		if not "{row_col}_colors".format(row_col = row_colors) in data.uns :
			useful.factor_colors(data, row_colors)

		row_colors = [data.uns["{row_col}_colors".format(row_col = row_colors)][gene] for gene in data.obs[row_colors]]

	if not "linkage" in data.uns :

		linkage = stat.calculate_linkage(data.X, data.obs.index.values)
		data.uns["linkage"] = linkage
	else :
		linkage = data.uns["linkage"]
	plot_df = pd.DataFrame(data.X, index = data.obs.index.values)
	plot = sns.clustermap(plot_df, row_linkage=linkage, col_linkage=linkage, row_colors = row_colors, **kwargs)
	return(plot)

def plot_kpca(
	data, 
	by_factor: str = None,
	label_factors: list = [],
	x_label: str = None,
	y_label: str = None,
	fig_height: float = 10,
	fig_width: float = 10,
	title: str = None,
	ax = None,
	**kwargs) :
	
	# Check the inputs
	if by_factor :
		if by_factor not in data.obs.columns :
			raise ValueError("{by_fac} not in .obs".format(by_fac = by_factor))

	if len(label_factors) > 0:
		if any([factor not in set(data.obs[by_factor]) for factor in label_factors]) :
			raise ValueError("{lab} not in {by_fac}".format(lab = label_factors, by_fac = by_factor))

	if not "k_pca" in data.obsm :
		print("No k-PCA found, calculating it")
		data.kernel_pca()

	pca = data.obsm["k_pca"]
	cols = None
	labels = None
	# Get colors for factor
	if by_factor :
		if not "{by_factor}_colors".format(by_factor = by_factor) in data.uns :
			useful.factor_colors(data, by_factor)
			col_dict = data.uns["{by_factor}_colors".format(by_factor = by_factor)]
		else :
			col_dict = data.uns["{by_factor}_colors".format(by_factor = by_factor)]

		if len(label_factors) > 0 :
			cols = [col_dict[factor] if factor in label_factors else "#D0D0D0" for factor in data.obs[by_factor]]
		else :
			cols = [col_dict[factor] for factor in data.obs[by_factor]]
		labels = [factor if factor in label_factors else None for factor in data.obs[by_factor]]
		
	plot = _scatter(x = pca[:,0], y = pca[:,1], c = cols, labels = labels,
	x_label = x_label, y_label = y_label, fig_height = fig_height, fig_width = fig_width, title = title, ax = ax, **kwargs)
	return(plot)

def plot_umap(
	data,
	by_factor: str = None,
	label_factors: list = [],
	x_label: str = None,
	y_label: str = None,
	fig_height: float = 10,
	fig_width: float = 10,
	title: str = None,
	**kwargs) :

	# Check the inputs
	if by_factor :
		if by_factor not in data.obs.columns :
			raise ValueError("{by_fac} not in .obs".format(by_fac = by_factor))

	# If True, label all factors
	if label_factors == True :
		label_factors = set(data.obs[by_factor])
	if len(label_factors) > 0:
		if any([factor not in set(data.obs[by_factor]) for factor in label_factors]) :
			raise ValueError("{lab} not in {by_fac}".format(lab = label_factors, by_fac = by_factor))

	if not "X_umap" in data.obsm :
		print("No UMAP found, calculating with 50 PCs")
		data.calc_umap()

	umap = data.obsm["X_umap"]
	
	cols = None
	labels = None

	if by_factor :
		if not "{by_factor}_colors".format(by_factor = by_factor) in data.uns :
			useful.factor_colors(data, by_factor)
			col_dict = data.uns["{by_factor}_colors".format(by_factor = by_factor)]
		else :
			col_dict = data.uns["{by_factor}_colors".format(by_factor = by_factor)]

		if len(label_factors) > 0 :
			cols = [col_dict[factor] if factor in label_factors else "#D0D0D0" for factor in data.obs[by_factor]]
		else :
			cols = [col_dict[factor] for factor in data.obs[by_factor]]
		labels = [factor if factor in label_factors else None for factor in data.obs[by_factor]]

	plot = _scatter(x = umap[:,0], y = umap[:,1], c = cols, labels = labels,
	x_label = x_label, y_label = y_label, fig_height = fig_height, fig_width = fig_width, title = title, **kwargs)
	
	return(plot)


def _scatter(
	x,
	y,
	labels = None,
	point_size: float = 5,
	s = None,
	c = None,
	x_label: str = None,
	y_label: str = None,
	fig_height: float = 10,
	fig_width: float = 10,
	title: str = None,
	ax = None,
	**kwargs) :

	# Need to group things by unique labels
	plot_dict = {}
	if labels :
		for key in set(labels) :
			x_vals = [val for val, label in zip(x, labels) if label == key]
			y_vals = [val for val, label in zip(y, labels) if label == key]
			c_vals = [val for val, label in zip(c, labels) if label == key]
			plot_dict[key] = {"x" : x_vals, "y" : y_vals, "c" : c_vals}
	else :
		plot_dict[None] = {"x" : x, "y" : y, "c" : c}

	if not ax :
		fig, ax = plt.subplots(1, 1, figsize = (fig_height, fig_width))
	for key, val in plot_dict.items() :
		ax.scatter(val["x"], val["y"], s = s, c = val["c"], label = key, **kwargs)

	if labels and c :
		ax.legend()
	ax.set_xlabel(x_label, fontsize=18)
	ax.set_ylabel(y_label, fontsize=16)
	ax.set_title(title, fontsize = 20)

	return(ax)

def plot_kpca_var(data, nPC = 100, **kwargs) :
	explained_variance = np.var(data.obsm["k_pca"], axis=0)
	explained_variance_ratio = explained_variance / np.sum(explained_variance) * 100
	explained_variance_ratio = explained_variance_ratio[0:nPC - 1]

	plot = plt.bar(x = range(len(explained_variance_ratio)), height = explained_variance_ratio, **kwargs)
	plot = plt.title("kpca explained variance")
	plot = plt.ylabel("percent of variance")
	plot = plt.xlabel("PC")

	return(plot)

def cluster_logo_plot(adata, obs_col, obs_val, lengths = "all") :
	length_args = ["all", "dominant"]
	if lengths not in length_args :
		raise ValueError("length argument must be one of %s" % length_args)

	# Lets take an example cluster
	logo_clust = adata.obs[adata.obs[obs_col] == obs_val][["TRA_cdr3", "TRB_cdr3", "TRA_cdr3_length", "TRB_cdr3_length"]]

	# Figure out the dominant lengths of the clusters
	if lengths == "dominant" :
		num_alpha_lengths = 1
		num_beta_lengths = 1
	else :
		num_beta_lengths = len(set(logo_clust["TRB_cdr3_length"]))
		num_alpha_lengths = len(set(logo_clust["TRA_cdr3_length"]))

	figRows = max([num_beta_lengths, num_alpha_lengths])

	# NEED TO FIGURE OUT HOW TO MAKE A GOOD FIGSIZE CALCULATION
	fig, ax = plt.subplots(nrows = figRows, ncols = 2, figsize = (10 * figRows, 3 * figRows))

	for num, seqtype in enumerate(["TRA", "TRB"]) :
		seq_df = logo_clust[["{seqtype}_cdr3".format(seqtype = seqtype), "{seqtype}_cdr3_length".format(seqtype = seqtype)]]

		if lengths == "dominant" :
			chain_lengths = [seq_df["{seqtype}_cdr3_length".format(seqtype = seqtype)].value_counts().idxmax()]
		else :
			chain_lengths = sorted(set(seq_df["{seqtype}_cdr3_length".format(seqtype = seqtype)]))
		for row, seqlen in enumerate(chain_lengths) :
			# Get the seqs
			logo_seqs = seq_df[seq_df["{seqtype}_cdr3_length".format(seqtype = seqtype)] == seqlen]["{seqtype}_cdr3".format(seqtype = seqtype)]
			# Concatenate and determine all used AA in seqs
			unique_AA = set("".join(seq for seq in logo_seqs))
			# Probability matrix
			prob_df = pd.DataFrame(index = range(seqlen), columns = unique_AA)

			for indx in range(len(logo_seqs[0])) :
				# Get the letter at the position for each seq
				AAs = [seq[indx] for seq in logo_seqs]

				# Calculate probabilities
				prob_dict = dict(Counter(AAs))
				for key, val in prob_dict.items() :
					prob = val / len(AAs)
					prob_df.loc[indx, key] = prob
			prob_df = prob_df.fillna(0)
			prob_df.sum(axis = 1)

			if figRows == 1 :
				logomaker.Logo(prob_df, ax = ax[num],
								width=.8,
								vpad=.05,
								fade_probabilities=True,
								stack_order='small_on_top',
								color_scheme='dodgerblue',
								font_name='Rosewood Std'
							   )
				ax[num].set_title("Number of seqs: {seqlen}".format(seqlen = len(logo_seqs)), {"fontsize" : 10, "fontweight" : "bold"})

				# Add additional title
				# Get the center of the plot
				center = seqlen / 1.75
				height = 1.1 + (figRows / 15)
				ax[num].text(center, height, "{seqtype} CDR3".format(seqtype = seqtype),{"fontsize" : 15, "fontweight" : "bold"} ,
								  horizontalalignment = "right")
				continue
			else :
				logomaker.Logo(prob_df, ax = ax[row, num],
								width=.8,
								vpad=.05,
								fade_probabilities=True,
								stack_order='small_on_top',
								color_scheme='dodgerblue',
								font_name='Rosewood Std'
							   )
				ax[row, num].set_title("Number of seqs: {seqlen}".format(seqlen = len(logo_seqs)), {"fontsize" : 10, "fontweight" : "bold"})
			# If the first of either alpha or beta, add additional title
				if row == 0 :
					# Get the center of the plot
					center = (seqlen + .75) / 2
					height = 1 + (figRows / 15)
					ax[row, num].text(center, height, "{seqtype} CDR3".format(seqtype = seqtype),{"fontsize" : 15, "fontweight" : "bold"} ,
									  horizontalalignment = "right")
	fig.tight_layout()

	# Determine which chain has more seq lengths
	if num_beta_lengths > num_alpha_lengths :
		# get rid of excess alpha plots
		for noLogo in range(num_alpha_lengths, num_beta_lengths) :
			ax[noLogo, 0].axis("off")
	elif num_beta_lengths < num_alpha_lengths :
		# Get rid of excess beta plots
		for noLogo in range(num_beta_lengths, num_alpha_lengths) :
			ax[noLogo, 1].axis("off")

	return(fig)