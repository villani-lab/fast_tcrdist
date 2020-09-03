import pandas as pd 
import numpy as np
import h5py
import time
import anndata
import scanpy as sc
def aggregate_cellranger_tcr(
	sample_sheet,
	save_location,
	chains = ["TRA", "TRB"],
	chain_vars = ["cdr3", "v_gene"]) :
	
	# Read in the sample info
	samp_info = pd.read_csv(sample_sheet)

	# Make sure columns are correct
	if not all([col in samp_info.columns for col in ["h5_location", "tcr_location"]]) :
		raise(ValueError( "sample sheet needs 2 columns : h5_location and tcr_location"))

	# Make a dict of chain/dtype combos
	# chain_dict = {key: chain_vars for key in chains}

	samp_dict = dict(zip(samp_info["h5_location"], samp_info["tcr_location"]))

	anndatas = []
	for gex, tcr in samp_dict.items() :
		gex_data = sc.read_10x_h5(gex)
		gex_data.var_names_make_unique()
		tcr_data = pd.read_csv(tcr)

		tcr_data = tcr_data[(tcr_data["productive"] == "True") & (tcr_data["high_confidence"] == True)]

		for dtype in chain_vars :
			for chain in chains :
				chain_data = tcr_data[tcr_data["chain"] == chain]
				info_dict = dict(zip(chain_data["barcode"], chain_data[dtype]))

				gex_data.obs["_".join([chain, dtype])] = [info_dict[code] if code in info_dict.keys() else "" for code in gex_data.obs.index]
		anndatas.append(gex_data)
	if len(anndatas) > 1 :
		all_data = anndatas[0].concatenate(anndatas[1:])
	else :
		all_data = anndatas[0]

	# Save the data
	all_data.write_h5ad(save_location)