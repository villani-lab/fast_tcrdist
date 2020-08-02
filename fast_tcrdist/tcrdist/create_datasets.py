import pandas as pd
import json
import os
from .. import plot
from .cython import seq_dist
from . import stats
from anndata import AnnData
import numpy as np
import time

class tcrdist_dataset(AnnData) :
    def __init__(self, seq_data, organism, only_betas = False) :
        super(tcrdist_dataset, self).__init__(shape = (seq_data.shape[0], seq_data.shape[0]))
        self.obs = seq_data
        self.obs.columns = self.obs.columns.str.replace(" ", "_")
        self.obs.index = self.obs.index.map(str)
        self.uns["species"] = organism
        self.uns["only_betas"] = only_betas

        if organism not in ["human", "mouse"] :
            raise(ValueError("organism must be either human or mouse"))
        with open(os.path.join(os.path.dirname(__file__), "db", "{organism}_v_gene_dists.json".format(organism = organism)), "r") as f:
            self.uns["v_gene_dists"] = json.load(f)
        for chain in ["TRA_cdr3", "TRB_cdr3"] :
            try :
                self.obs["_".join([chain, "trimmed"])] = self.trim_cdr3s(seq_data[chain])
            except :
                continue

        # Get information about CDR3 length
        for chain in ["TRA_cdr3", "TRB_cdr3"] :
            try :
                self.obs["_".join([chain, "length"])] = self.obs.apply(lambda df : len(df[chain]), axis = 1)
            except :
                continue

    def trim_cdr3s(self, seqs, front_trim = 3, back_trim = 2) :
        '''tcrdist calculation used a trimmed CDR3:
         "trimmed CDR3’ defined as starting with the 3rd position after the C104 
         and terminating with the 2nd position before F118 (5 fewer residues than the full CDR3"
        Dash et al. 2017.
         '''
        trimmed_seqs = [seq[front_trim:-back_trim] for seq in seqs]
        return(trimmed_seqs)

    def calc_distance_legacy(self, seq1_info, seq2_info) :
        
        # Calculate the distance of the CDR3s
        aligned_A = seq_dist.align_cdrs(seq1_info["TRA_cdr3_trimmed"], seq2_info["TRA_cdr3_trimmed"])
        cdr3_dist_A = seq_dist.blosum_sequence_distance(aligned_A[0], aligned_A[1], gap_penalty = 12, cdr3 = True)

        aligned_B = seq_dist.align_cdrs(seq1_info["TRB_cdr3_trimmed"], seq2_info["TRB_cdr3_trimmed"])
        cdr3_dist_B = seq_dist.blosum_sequence_distance(aligned_B[0], aligned_B[1], gap_penalty = 12, cdr3 = True)

        # Get the distances of the V genes
        v_dist_A = self.uns["v_gene_dists"]["_".join([seq1_info["TRA_v_gene"], seq2_info["TRA_v_gene"]])]
        v_dist_B = self.uns["v_gene_dists"]["_".join([seq1_info["TRB_v_gene"], seq2_info["TRB_v_gene"]])]
        tot_dist = sum([cdr3_dist_A, cdr3_dist_B, v_dist_A, v_dist_B])

        return(tot_dist)

    def create_tcrdist_matrix_legacy(self, barcodes) :
        import time
        start = time.time()
        mtx = np.zeros((len(barcodes), len(barcodes)))

        for code1 in range(len(barcodes)) :
            for code2 in range(code1 + 1, len(barcodes)) :
                data1 = self.collect_info(barcodes[code1])
                data2 = self.collect_info(barcodes[code2])

                dist = self.calc_distance_legacy(data1, data2)
                mtx[code1, code2] = dist
                mtx[code2, code1] = dist
        end = time.time()
        print("time for mtx: %s seconds"%(end - start))
        self.X = mtx
        return(mtx)
    
    def create_tcrdist_matrix_nw(self) :
            start = time.time()
            mtx = seq_dist.create_tcrdist_matrix(self.obs, v_scores = self.uns["v_gene_dists"], only_betas = self.uns["only_betas"])
            end = time.time()
            print("time for mtx: %s seconds"%(end - start))
            print(mtx)
            self.X = mtx

    def collect_info(self, barcode) :
        data = self.obs.loc[barcode]
        info = {
        "TRA_cdr3_trimmed" : data["TRA_cdr3_trimmed"],
        "TRB_cdr3_trimmed" : data["TRB_cdr3_trimmed"],
        "TRA_v_gene" : data["TRA_v_gene"],
        "TRB_v_gene" : data["TRB_v_gene"]
        }
        return(info)

    def kernel_pca(self, kernel = "precomputed") :
        pca = stats.calc_kernel_PCA(self.X, kernel = kernel)
        self.obsm["k_pca"] = pca

    def calc_umap(self, return_coords = False) :
        import umap
        # Make sure kernel-PCA has been calculated
        # if not "k_pca" in self.obsm :
        #     raise ValueError("no PCA in object.  Please run kernel_pca first!")
        
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(self.X)

        self.obsm["X_umap"] = embedding

        if return_coords :
            return(embedding)

    def hdbscan_cluster(self, return_results = False, coord = "umap", n_coord = 2, **kwargs) :
        embedding = "_".join(["X", coord])
        if not embedding in self.obsm :
            raise ValueError("No {coord} in object.  Please run {coord} first!".format(coord = coord))
            
        clust = stats.hdbscan(self.obsm[embedding][:,0:n_coord], **kwargs)
        # Need to add one because pandas gets confused with categorical data < 1, make -1 "unassigned"
        clust = [str(num + 1) if num != -1 else "unassigned" for num in clust]
        self.obs["hdbscan_clust"] = clust

        if return_results :
            return(clust)

def tcrdist_dataset_from_anndata(adata, organism, obs = []) :
    '''Function to initialize an object of 'tcrdist_dataset' from an anndata that has been processed via pegasus'''

    adata = adata[(adata.obs["TRB_cdr3"] != "") & (adata.obs["TRA_cdr3"] != "") & (adata.obs["TRA_v_gene"] != "") & (adata.obs["TRB_v_gene"] != "")]
    
    # Get the needed columns
    need_cols = ["TRA_cdr3", "TRB_cdr3", "TRA_v_gene", "TRB_v_gene"] + obs

    needed_data = adata.obs[need_cols]
    # Removes categorical datatype if it exists
    needed_data = pd.DataFrame(needed_data, dtype = "object")
    dset = tcrdist_dataset(seq_data = needed_data, organism = organism)

    return(dset)

def tcrdist_dataset_from_vdjdb(tsv_path, organism, only_betas = False, cell_id = "complex.id", chain_col = "Gene", cdr3_col = "CDR3", v_gene_col = "V",
    other_cols = []) :
    '''Function to initialize a tcrdist object from data downloaded from vdjdb (vdjdb.cdr3.net)'''
    data = pd.read_csv(tsv_path, sep = "\t")

    # Get just the data we need
    data = data[[cell_id, chain_col, cdr3_col, v_gene_col] + other_cols]
    # Adjust the formatting of the V genes
    data["V"] = [str(gene).replace("*", "-") for gene in data["V"]]
    # remove additional allele information(not needed)
    data["V"] = [gene[:-3] for gene in data["V"]]
    # Get rid of the slashes  in some of the names
    data["V"] = [gene.replace("/", "") for gene in data["V"]]

    if only_betas :
        needed_data = data[data["Gene"] == "TRB"]
        needed_data = needed_data.rename(columns = {"CDR3" : "TRB_cdr3", "V" : "TRB_v_gene"})
        needed_data = needed_data.drop(["Gene", "complex.id"], axis = 1)

    else :
    # Seperate the alpha and beta data
        df_dict = {}
        for chain in ["TRA", "TRB"] :
            df = data[data["Gene"] == chain]
            if chain == "TRA" :
                df = df.drop(other_cols, axis = 1)
            df = df.rename(columns = {"CDR3" : "{chain}_cdr3".format(chain = chain), "V" : "{chain}_v_gene".format(chain = chain)})
            df = df.drop(["Gene"], axis = 1)
            df = df.set_index("complex.id")

            df_dict[chain] = df
        needed_data = df_dict["TRA"].join(df_dict["TRB"], on = "complex.id")

        # Remove the TCRs with missing information
        needed_data = needed_data[(needed_data["TRB_cdr3"] != "") &
        (needed_data["TRA_cdr3"] != "") & 
        (needed_data["TRA_v_gene"] != "") & 
        (needed_data["TRB_v_gene"] != "")]
        needed_data = needed_data.dropna()

    dset = tcrdist_dataset(seq_data = needed_data, organism = organism)
    return(dset)

def tcrdist_dataset_from_csv(csv_path, organism) :
    # needs to be done
    return(0)

def tcrdist_dataset_from_df(df, organism, only_betas = False, from_adaptive = False, cdr3a_col = "TRA_crd3", cdr3b_col = "TRB_cdr3", Va_col = "TRA_v_gene", Vb_col = "TRB_v_gene") :
    if only_betas :
        needed_data = df
        needed_data = needed_data.rename(columns = {cdr3b_col : "TRB_cdr3", Vb_col : "TRB_v_gene"})
    else :
        needed_data = df
        needed_data = needed_data.rename(columns = {cdr3a_col : "TRA_cdr3", cdr3b_col : "TRB_cdr3", Va_col: "TRA_v_gene", Vb_col : "TRB_v_gene"})

    if from_adaptive :
        needed_data["TRB_v_gene"] = [gene if gene != "TRBV9-1" else "TRBV9" for gene in needed_data["TRB_v_gene"]]
        needed_data["TRB_v_gene"] = [gene if gene != "TRBV26-1" else "TRBV26" for gene in needed_data["TRB_v_gene"]]
        needed_data["TRB_v_gene"] = [gene if gene != "TRBV1-1" else "TRBV1" for gene in needed_data["TRB_v_gene"]]

    dset = tcrdist_dataset(seq_data = needed_data, organism = organism, only_betas = only_betas)
    return(dset)