def random_colors (n_cols) :
	import random
	num = lambda: random.randint(0,255)

	cols = set()
	for i in range(n_cols) :
		cols.add("#%02X%02X%02X" % (num(),num(),num()))
	return(cols)

def create_color_dict (categories) :
	colors = random_colors(len(categories))

	color_dict = dict(zip(categories, colors))

	return(color_dict)

def factor_colors(data, factor):
	
	# Get the unique genes
	unique_factors = set(data.obs[factor])
	cdict = create_color_dict(unique_factors)

	data.uns["{factor}_colors".format(factor = factor)] = cdict