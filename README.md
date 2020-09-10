# fast_tcrdist

fast_tcrdist is an optimized version of the TCRDist algorithm published by Dash et al. Nature (2017): doi:10.1038/nature22383

To enhance the original implementation of TCRDist, fast_tcrdist uses the Needleman-Wunsch alignment algorithm to align TCR sequences
and creates the TCRDist matrix via cython.  To integrate well with other common single cell analysis tools, fast_tcrdist utilizes
the Anndata data structure to store the TCRDist matrix and associated metadata.  Currently, this has been tested on TCR/gene expression
output from Cellranger (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), but 
future releases will aim to allow for other file formats.

In addition to running the TCRDist algorithm, fast_tcrdist allows you to aggregate TCR info with single-cell gene expression into a single
anndata object to allow for integrated downstream analyses.

**Outside files**

The BLOSUM62 matrix used for alignments came from NCBI:https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM62

CDR amino acid info was taken from the TCRDist database file "alphabeta_db.tsv" (https://www.dropbox.com/s/kivfp27gbz2m2st/tcrdist_extras_v2.tgz)
and reformated into .json format (mouse_CDRs_for_10X.json and human_CDRs_for_10X.json)