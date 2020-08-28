import os
from . import blosum_mtx
cdef dict adj_blosum = blosum_mtx.adj_blosum
cdef dict blosum = blosum_mtx.blosum

import numpy as np
cimport numpy as np
from .cnwalign import global_align
matrix = os.path.join(os.path.dirname(os.path.abspath(__file__)), "BLOSUM62")

cpdef tuple align_cdrs(a, b):
	# If 2 CDR3s same length, directly align
	if len(a) == len(b):
		return (bytes(a[:], "utf-8"), bytes(b[:], "utf-8"))
	
	# Else assign the sequences to s0, s1 where s0 is the shorter one
	# cdef char* s0
	# cdef char* s1
	if len(a)<len(b): 
		s0,s1 = a,b
	else:
		s0,s1 = b,a
	
	# Determine difference in length
	cdef int lendiff = len(s1)-len(s0)
	
	cdef int best_score=-1000
	cdef int best_gappos=0 # in case len(s0) == 1

	# the gap comes after s0[gappos]
	# Iterate across all possible gaps in the sequences
	cdef int gappos
	cdef int score
	for gappos in range(len(s0)-1):
		score=0
		# Iterating across the CDR3s, aligning them to the left
		for i in range(gappos+1):
			# Determine the blosum score of the two AA, add to score
			score += blosum[s0[i] + s1[i]]
		# Iterating across the CDR3s,determining the scores when CDR3s are different lengths
		# using lendiff ensures you make it to the end of the second sequqence
		# Through all iterations, you will try all possible combinations
		for i in range(gappos+1,len(s0)):
			score += blosum[s0[i] + s1[i+lendiff]]
		# Determine if the score obtained is the best of all iterations
		if score>best_score:
			# If best score, take it
			best_score = score
			# Also take location of possible gap
			best_gappos = gappos
	# Adjust the shorter of the 2 sequences to now include the gap 
	# Gap represented by gap_character(currently ".")
	temp = s0[:best_gappos+1] + "-" * lendiff + s0[best_gappos+1:] 
	s0 = temp
	
	# Make sure length of the sequences are now the same
	assert len(s0) == len(s1)

	# Return the new sequences
	if len(a)<len(b): 
		return (bytes(s0, "utf-8"), bytes(s1, "utf-8"))
	else:
		return (bytes(s1, "utf-8"), bytes(s0, "utf-8"))

cpdef int blosum_character_distance(char* a, char* b, int gap_penalty, cdr3 = False):
	gap_char = b"-"
	if a== gap_char and b == gap_char :
		return 0
	elif a == gap_char or b == gap_char :
		return gap_penalty
	else:
		if cdr3 :
			return adj_blosum[a+b] * 3
		else :
			return adj_blosum[a+b]

cpdef int blosum_sequence_distance(char* aseq, char* bseq, int gap_penalty, cdr3 = False):
	assert len(aseq) == len(bseq)
	cdef int dist = 0
	cdef Py_ssize_t num
	cdef Py_ssize_t indx = len(aseq)
	for num in range(indx) :
		dist += blosum_character_distance(aseq[num: num + 1], bseq[num:num + 1], gap_penalty, cdr3 = cdr3)
	return(dist)

cpdef int nw_align(char* a, char* b) :
	cdef int score
	if len(a) == len(b) :
		score = blosum_sequence_distance(a, b, gap_penalty = 12, cdr3 = True)
		return(score)

	cdef tuple align = global_align(a, b, matrix = matrix, gap_open = -20, gap_extend = -20)
	
	score = blosum_sequence_distance(align[0], align[1], gap_penalty = 12, cdr3 = True)
	return(score)


cpdef np.ndarray create_tcrdist_matrix(dataframe, v_scores, only_betas = False) :
	import itertools
	# Create V-gene matrices
	cdef Py_ssize_t mtx_dim = dataframe.shape[0]
	cdef tuple triu = np.triu_indices(mtx_dim) # Find upper right indices of a triangular nxn matrix
	cdef tuple tril = np.tril_indices(mtx_dim, -1) # Find lower left indices of a triangular nxn matrix

	cdef np.ndarray mtxb = np.zeros((mtx_dim, mtx_dim), dtype = int)
	vb_combinations = itertools.combinations_with_replacement(dataframe["TRB_v_gene"], 2)
	vb_scores = [v_scores["_".join([seq1, seq2])] for seq1, seq2 in vb_combinations]
	mtxb[triu] = vb_scores # Assign list values to upper right matrix

	cdef list va_scores
	cdef np.ndarray mtxa = np.zeros((mtx_dim, mtx_dim), dtype = int)
	if not only_betas :
		va_combinations = itertools.combinations_with_replacement(dataframe["TRA_v_gene"], 2)
		va_scores = [v_scores["_".join([seq1, seq2])] for seq1, seq2 in va_combinations]

		mtxa[triu] = va_scores 


	# Make the CDR3 matrix
	cdef np.ndarray mtx_cdr3 = np.zeros((mtx_dim, mtx_dim), dtype = int)

	# Convert all CDR3s to bytes and put into an ndarray
	cdef np.ndarray betas = dataframe["TRB_cdr3_trimmed"].to_numpy(dtype = "<S36")
	cdr3b_combinations = itertools.combinations_with_replacement(betas, 2)
	cdef np.ndarray alignb = np.zeros(triu[0].shape[0], dtype = int)
	cdef int cntr = 0
	cdef char* seq_1
	cdef char* seq_2
	for seq_1, seq_2 in cdr3b_combinations :
		alignb[cntr] = nw_align(seq_1, seq_2)
		cntr +=1
	
	cdef np.ndarray aligna = np.zeros(triu[0].shape[0], dtype = int)
	cdef int cntr1
	cdef np.ndarray alphas
	if not only_betas :
		alphas = dataframe["TRA_cdr3_trimmed"].to_numpy(dtype = "<S36")
		cdr3a_combinations = itertools.combinations_with_replacement(alphas, 2)

		cntr1 = 0
		for seq_1, seq_2 in cdr3a_combinations :
			aligna[cntr1] = nw_align(seq_1, seq_2)
			cntr1 +=1

	cdef np.ndarray cdr3_scores = np.sum([alignb, aligna], axis = 0)

	mtx_cdr3[triu] = cdr3_scores

	cdef np.ndarray mtx = np.sum([mtxb, mtxa, mtx_cdr3], axis = 0)
	mtx[tril] = mtx.T[tril] # Make the matrix symetric

	return(mtx)

## Work in progress, implement NW in Cython to not rely on other packages ##
# cpdef int nw_align(char* x, char* y, gap = -20):
# 	nx = len(x)
# 	ny = len(y)
# 	# print(x)
# 	# Optimal score at each possible pair of characters.
# 	F = np.zeros((nx + 1, ny + 1))
# 	F[:,0] = np.linspace(0, -nx, nx + 1)
# 	F[0,:] = np.linspace(0, -ny, ny + 1)
# 	# Pointers to trace through an optimal aligment.
# 	P = np.zeros((nx + 1, ny + 1))
# 	P[:,0] = 3
# 	P[0,:] = 4
# 	# Temporary scores.
# 	t = np.zeros(3)
# 	for i in range(nx):
# 		for j in range(ny):
# 			# if x[i] == y[j]:
# 			#     t[0] = F[i,j] + blosum[x[i] + y[j]]
# 			# else:
# 			#     t[0] = F[i,j] + blosum[x[i] + y[j]]
# 			# Need to get the proper index

# 			t[0] = F[i, j] + blosum[x[i] + y[j]]

# 			t[1] = F[i,j+1] + gap
# 			t[2] = F[i+1,j] + gap
# 			tmax = np.max(t)
# 			F[i+1,j+1] = tmax
# 			if t[0] == tmax:
# 				P[i+1,j+1] += 2
# 			if t[1] == tmax:
# 				P[i+1,j+1] += 3
# 			if t[2] == tmax:
# 				P[i+1,j+1] += 4
# 	# Trace through an optimal alignment.
# 	i = nx
# 	j = ny
# 	rx = []
# 	ry = []
# 	while i > 0 or j > 0:
# 		if P[i,j] in [2, 5, 6, 9]:
# 			rx.append(x[i-1])
# 			ry.append(y[j-1])
# 			i -= 1
# 			j -= 1
# 		elif P[i,j] in [3, 5, 7, 9]:
# 			rx.append(x[i-1])
# 			ry.append('-')
# 			i -= 1
# 		elif P[i,j] in [4, 6, 7, 9]:
# 			rx.append('-')
# 			ry.append(y[j-1])
# 			j -= 1
# 	# Reverse the strings.
# 	rx = ''.join(rx)[::-1]
# 	ry = ''.join(ry)[::-1]

# 	rx = bytes(rx, "utf-8")
# 	ry = bytes(ry, "utf-8")

# 	# Make sure strings are same length
# 	assert len(rx) == len(ry)
# 	score = blosum_sequence_distance(rx, ry, gap_penalty = 12, cdr3 = True)
# 	# print(bytes(rx, "utf-8"))
# 	# print(bytes(ry, "utf-8"))
# 	# score = 0
# 	return(score)
