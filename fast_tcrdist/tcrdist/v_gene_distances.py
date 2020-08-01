'''Compute distances of all V genes and store for reference'''
import json
import distance

def compute_all_v_region_distances(organism):
    with open("db/{organism}_cdrs_for_10X.json".format(organism = organism), "r") as f:
        cdr_dict = json.load(f)

    # Get all unique V gene names
    keys = [list(cdr_dict[key].keys()) for key in cdr_dict.keys()]
    assert keys[0] == keys[1] == keys[2]
    keys = keys[0]

    # Collect all of the information for each V gene
    v_gene_dict = {}
    for v_gene in keys :
        # Get all of the cdrs for each V gene
        cdrs = [cdr_dict[cdr][v_gene] for cdr in cdr_dict.keys()]
        
        v_gene_dict[v_gene] = cdrs
    v_comp_dict = {}
    # Get distance of possible V gene matches
    for chain in ["A", "B"] :
        # Get V genes from a single chain
        seqs = {key : v_gene_dict[key] for key in v_gene_dict.keys() if chain in key}
        for v_1 in seqs.keys() :
            for v_2 in seqs.keys() :
                aligned = [distance.align_cdrs(seq1, seq2) for seq1, seq2 in zip(seqs[v_1], seqs[v_2])]
                dists = [distance.blosum_sequence_distance(seq1, seq2, gap_penalty = 4) for seq1, seq2 in aligned]
                total_dist = sum(dists)
                v_comp_dict["_".join([v_1, v_2])] = total_dist
    with open("db/{organism}_v_gene_dists.json".format(organism = organism), "w") as f :
            json.dump(v_comp_dict, f, indent = 4)
            f.close()

if __name__ == "__main__" :
    compute_all_v_region_distances("human")
    compute_all_v_region_distances("mouse")