import pandas as pd

def get_allowed_chrom_str():
	return ['chr%d'%(i+1) for i in range(22)]

def get_chrom_sizes(sizes_path):
	df = pd.read_csv(sizes_path, names=['chrom', 'size'], sep='\t')
	return dict(zip(df['chrom'],df['size']))
