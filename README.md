# fast_tcrdist

fast_tcrdist is an optimized version of the TCRDist algorithm published by Dash et al. Nature (2017): doi:10.1038/nature22383

To enhance the original implementation of TCRDist, fast_tcrdist uses the Needleman-Wunsch alignment algorithm to align TCR sequences
and creates the TCRDist matrix via cython.  To integrate well with other common single cell analysis tools, fast_tcrdist utilizes
the Anndata data structure to store the TCRDist matrix and associated metadata.