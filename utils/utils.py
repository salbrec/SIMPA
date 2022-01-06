import pandas as pd

def get_allowed_chrom_str(assembly):
	chroms = None
	if assembly in ['hg19', 'hg38']:
		chroms = ['chr%d'%(i+1) for i in range(22)]
	if assembly in ['mm10']:
		chroms = ['chr%d'%(i+1) for i in range(19)]
	return chroms

def get_chrom_sizes(sizes_path):
	df = pd.read_csv(sizes_path, names=['chrom', 'size'], sep='\t')
	return dict(zip(df['chrom'],df['size']))
