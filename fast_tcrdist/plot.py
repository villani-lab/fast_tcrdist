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
	row_colors : str = None,
	**kwargs) :
	
	if row_colors :
		if row_colors not in  data.obs.columns :
			raise ValueError("{row_col} not in .obs".format(row_col = row_colors))
		# Get the categories
		if not pd.api.types.is_categorical(data.obs[row_colors]) :
			data.obs[row_colors] = data.obs[row_colors].astype("category")
		categories = data.obs[row_colors].cat.categories

		if not "{row_col}_colors".format(row_col = row_colors) in data.uns :

			palette = useful.godsnot_102[:len(categories)]
			data.uns["{row_col}_colors".format(row_col = row_colors)] = palette

		col_dict = dict(zip(categories, data.uns["{row_col}_colors".format(row_col = row_colors)]))

		row_colors = [col_dict[gene] for gene in data.obs[row_colors]]

	if not "linkage" in data.uns :
		linkage = stat.calculate_linkage(data.X, data.obs.index.values)
		data.uns["linkage"] = linkage
	else :
		linkage = data.uns["linkage"]
	plot_df = pd.DataFrame(data.X, index = data.obs.index.values)
	plot = sns.clustermap(plot_df, row_linkage=linkage, col_linkage=linkage, row_colors = row_colors, **kwargs)
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