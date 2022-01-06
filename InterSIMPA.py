import pandas as pd
import numpy as np
from sys import argv, path, stdout
import sys
import os
import pickle
import argparse
from terminaltables import AsciiTable
from sklearn.ensemble import RandomForestClassifier

import utils.pp as pp
import utils.utils as utils

def print_nice_table(table):
	print_table = AsciiTable(table)
	print(print_table.table)

simpa_dir = './'
if argv[0].find('/') >= 0:
	simpa_dir = argv[0][: - argv[0][::-1].find('/')]

parser = argparse.ArgumentParser(description='InterSIMPA - INTERpretation for imputed probability of SIMPA')
parser.add_argument('--bed', '-b', type=str, required=True, help='Path to bed file with sparse single-cell input')
parser.add_argument('--targets', '-t', type=str, required=True, help='''Target(s) defining the specific reference experiments
					(ususally the one used in the scChIP). When multiple targets are provided, separate by "+"''')
parser.add_argument('--summit', '-s', type=str, required=True, help='Peak summit of target region to be investigated')
parser.add_argument('--outdir', '-o', type=str, default='./', help='Output directory. Default: "./"')
parser.add_argument('--genome', '-g', type=str, default='hg38', choices=['hg38'], help='Genome assembly')
parser.add_argument('--binsize', '-bs', type=str, default='5kb', choices=['5kb','50kb'], help='Size of the bins (genomic regions)')
parser.add_argument('--estimators', '-e', type=int, default=1000, help='Number of trees in Random Forest')
parser.add_argument('--importance', '-it', type=float, default=1.0, help='Threshold for the feature importance')


# parse and pre-process command line arguments
args = parser.parse_args()
bin_size_int = int(args.binsize.replace('bp', '').replace('kb', '000'))
ENCODE_dir = '%sdata/ENCODE/%s/%s/'%(simpa_dir, args.genome, args.binsize)
target_set = set(args.targets.split('+'))

# initialize variables used by all ranks
sc_bins = None
training_features = None
uniq_cand_labels = None
candidate = None
freq_map = None
ref_n_bins = None
bin_bed_map = None

# used to extract uniq bin columns
bin_index_map = None
key_index_map = None

# initialize variables needed to read in the sparse single-cell input
# that is converted into a set (or list) of bins
allowed_chroms = utils.get_allowed_chrom_str(args.genome)
chrom_sizes = utils.get_chrom_sizes('%sdata/chromosome_sizes/%s.tsv'%(simpa_dir, args.genome))
peaks = pp.get_peaks(args.bed, allowed_chroms, enrich_index=-1)
sc_bins, sc_bin_value, max_bin_ID, bin_bed_map = pp.bin_it(peaks, allowed_chroms,
														chrom_sizes, bin_size_int)
sc_bins = sorted(list(sc_bins))
print('\nGiven the sparse input there are %d genomic regions converted into %d bins of size %s'%(
	sum([len(p) for p in peaks.values()]), len(sc_bins), args.binsize))

# get reference experiments for given target(s)
metadata = pd.read_csv('%sdata/metadata_ENCODE.tsv'%(simpa_dir), sep='\t')
metadata = metadata.loc[metadata['assembly'] == args.genome]
metadata = metadata.loc[[True if target in target_set else False for target in metadata['target']]]
print('Number of available bulk reference experiments: %d (for %s)'%(metadata.shape[0],args.targets))

# read in the reference experiments, already converted into bins
# collect also all bins that are observed in at least one reference experiment
all_ref_bins = set()
ref_bins_map = {}
for accession in metadata['accession']:
	ref_experiment_bins = pickle.load(open(ENCODE_dir+accession, 'rb'))
	all_ref_bins = all_ref_bins | ref_experiment_bins
	ref_bins_map[accession] = ref_experiment_bins

# calculate the frequencies: how often is a particular bin present across all the
# reference experiments
#print('Number of bins with a signal in at least one bulk:', len(all_ref_bins))
freq_map = {}
for bid in all_ref_bins:
	freq = sum([True if bid in bin_set else False for bin_set in ref_bins_map.values()])
	freq = float(freq) / float(metadata.shape[0])
	freq_map[bid] = freq

# number of bins within a reference set
ref_n_bins = [len(ref_bins_map[acc]) for acc in metadata['accession']]

cand_chr = args.summit.split(':')[0]
cand_start = int(args.summit.split(':')[1])
cand_end = cand_start + 1

cand_peak = {cand_chr:[(cand_chr, cand_start, cand_end, -1)]}
cand_bins, cand_bin_value, max_bin_ID_2, bin_bed_map_2 = pp.bin_it(cand_peak, allowed_chroms,
														chrom_sizes, bin_size_int)

cand_bin = list(cand_bins)[0]
print('\nPeak summit for given target region:', cand_chr, cand_start, cand_end)
print('Bin ID for given target region:', cand_bin)

bin_is_present = cand_bin in sc_bins
if bin_is_present:
	sc_bins.remove(cand_bin)

# depending on the bins in the sparse input, create the training features matrix
training_features = np.zeros((metadata.shape[0], len(sc_bins)), dtype=np.int)
for index, accession in enumerate(metadata['accession']):
	training_features[index] = np.array([1 if bid in ref_bins_map[accession] else 0 for bid in sc_bins])
training_features = np.array(training_features)
print('Shape of matrix for training features:', training_features.shape)
print('\n##### Pre-Processing is done ... #####\n')

class_vector = np.array([1 if cand_bin in ref_bins_map[acc] else 0 for acc in metadata['accession']])
ref_frequency = np.mean(class_vector)

print('Bin is present in given cell:', 'YES' if bin_is_present else 'NO')
print('Reference Frequency of the given bin: %.3f'%(ref_frequency))

if ref_frequency > 0.9999999:
	print('\nAs the reference frequency is 100%, it is not possible to train and interpret a model')
	sys.exit()

if ref_frequency < 0.0000001:
	print('\nAs the reference frequency is 0%, it is not possible to train and interpret a model')
	sys.exit()

# initialize classification algorithm with Random Forest
clf = RandomForestClassifier(n_estimators=args.estimators, random_state=42)

clf.fit(training_features, class_vector)
artificial_inst = np.array([training_features.shape[1] * [1]])
prob = clf.predict_proba(artificial_inst)[0][1]

print('The imputed probability: %.3f\n'%(prob))

sc_bin_importance = list(zip(sc_bins, clf.feature_importances_ * 100.0))
sc_bin_importance = sorted(sc_bin_importance, key=lambda x: x[1], reverse=True)

gene_info = pd.read_csv('%sdata/genes/THE_GENES.tsv'%(simpa_dir), sep='\t', low_memory=False)
gene_descr = dict(zip(gene_info['geneID'], gene_info['description']))
gene_symbol = dict(zip(gene_info['geneID'], gene_info['symbol']))
gene_orient = dict(zip(gene_info['geneID'], gene_info['orientation']))
gene_chrom = dict(zip(gene_info['geneID'], gene_info['chromosome']))
gene_start = dict(zip(gene_info['geneID'], gene_info['start_pos']))
gene_end = dict(zip(gene_info['geneID'], gene_info['end_pos']))

genes_pos = dict( (chrom, []) for chrom in allowed_chroms )
for gene_id in gene_info['geneID']:
	if not gene_chrom[gene_id] in allowed_chroms:
		continue
	if np.isnan(gene_start[gene_id]) or np.isnan(gene_end[gene_id]):
		continue

	tss = gene_start[gene_id] if gene_orient[gene_id] == 'plus' else gene_end[gene_id]
	genes_pos[gene_chrom[gene_id]].append( (gene_id, int(gene_start[gene_id]), int(gene_end[gene_id]), tss) )

print('\nFeatures, represented by bins, most important for the model (importance > 1%)')
print('additional annotations regarding the next gene (if bin overlaps gene-body, Distance==0)\n')

table = [['Bin ID','Importance','Genomic Region','Next Gene (NG)','Dist TSS',
		  'Genebody','Dist TTS','NG Orientation','NG Description']]
for bin_id, importance in sc_bin_importance:
	if importance < args.importance:
		break
	bin_chrom = bin_bed_map[bin_id][0]
	bin_start = bin_bed_map[bin_id][1]
	bin_end   = bin_bed_map[bin_id][2]
	bin_mid   = int(np.mean( [bin_start, bin_end] ))

	next_gene = None
	gene_body = False

	# first search for gene bodies that overlap the bin
	for gene_id, start, end, tss in genes_pos[bin_chrom]:
		if end > bin_start:
			if start < bin_end:
				next_gene = (gene_id, 0)
				gene_body = True
				break

	# if there is no overlapping gene body:
	# annotate region by closest TSS
	if next_gene == None:
		genes_dist = [ (gene_id, min( abs(tss - bin_start), abs(tss - bin_end) )) for gene_id, start, end, tss in genes_pos[bin_chrom] ]
		next_gene = sorted(genes_dist, key=lambda x: x[1])[0]

	next_gene_id, next_gene_dist = next_gene

	# get further information about the distance to the gene
	tss = gene_start[next_gene_id] if gene_orient[next_gene_id] == 'plus' else gene_end[next_gene_id]
	tts = gene_end[next_gene_id] if gene_orient[next_gene_id] == 'plus' else gene_start[next_gene_id]
	dist_tss = int(bin_mid - tss)
	dist_tts = int(bin_mid - tts)
	if gene_orient.get(next_gene_id,'???') == 'minus':
		dist_tss *= -1
		dist_tts *= -1

	row = [bin_id, '%.1f%s'%(importance,'%'), '%s:%d-%d'%(bin_bed_map[bin_id])]
	row += [gene_symbol[next_gene_id], str(dist_tss), 'overlap' if gene_body else 'no',
		 str(dist_tts), gene_orient.get(next_gene_id,'???'), gene_descr.get(next_gene_id,'???')]
	table.append(row)

print_nice_table(table)
print('')



















