import pandas as pd
import numpy as np
from sys import argv, path, stdout
import os
import pickle
import argparse
from sklearn.ensemble import RandomForestClassifier

# initializing MPI
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import utils.pp as pp
import utils.utils as utils

spartan_dir = './'
if argv[0].find('/') >= 0:
	spartan_dir = argv[0][: - argv[0][::-1].find('/')]

parser = argparse.ArgumentParser(description='SPARTAN: SPARse peaks impuTAtioN')
parser.add_argument('--bed', '-b', type=str, required=True, help='Path to bed file with sparse single-cell input')
parser.add_argument('--targets', '-t', type=str, required=True, help='''Target(s) defining the specific reference experiments 
					(ususally the one used in the scChIP). When multiple targets are provided, separate by "+"''')
parser.add_argument('--outdir', '-o', type=str, default='./', help='Output directory. Default: "./"')
parser.add_argument('--genome', '-g', type=str, default='hg38', choices=['hg38'], help='Genome assembly')
parser.add_argument('--binsize', '-bs', type=str, default='5kb', choices=['5kb','50kb'], help='Size of the bins (genomic regions)')
parser.add_argument('--estimators', '-e', type=int, default=100, help='Number of trees in Random Forest')
parser.add_argument('--simulate', action='store_true', help='Impute only 100 bins for testing the software')

# parse and pre-process command line arguments
args = parser.parse_args()
bin_size_int = int(args.binsize.replace('bp', '').replace('kb', '000'))
ENCODE_dir = '%sdata/ENCODE/%s/%s/'%(spartan_dir, args.genome, args.binsize)
target_set = set(args.targets.split('+'))

# initialize variables used by all ranks
sc_bins = None
imp_prob_all = None
dimensions = None
training_features = None
uniq_cand_labels = None
candidates = None
freq_map = None
ref_n_bins = None
bin_bed_map = None

# used to extract uniq bin columns
bin_index_map = None
key_index_map = None

if rank == 0:
	# initialize variables needed to read in the sparse single-cell input
	# that is converted into a set (or list) of bins
	allowed_chroms = utils.get_allowed_chrom_str()
	chrom_sizes = utils.get_chrom_sizes('%sdata/chromosome_sizes/%s/sizes.tsv'%(spartan_dir, args.genome))
	peaks = pp.get_peaks(args.bed, allowed_chroms, enrich_index=-1)
	sc_bins, sc_bin_value, max_bin_ID, bin_bed_map = pp.bin_it(peaks, allowed_chroms, 
															chrom_sizes, bin_size_int)
	sc_bins = sorted(list(sc_bins))
	print('\nGiven the sparse input there are %d genomic regions converted into %d bins of size %s'%(
		sum([len(p) for p in peaks.values()]), len(sc_bins), args.binsize))
	
	# get reference experiments for given target(s)
	metadata = pd.read_csv('%sdata/metadata_ENCODE.tsv'%(spartan_dir), sep='\t')
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
	print('Number of bins with a signal in at least one bulk:', len(all_ref_bins))
	freq_map = {}
	for bid in all_ref_bins:
		freq = sum([True if bid in bin_set else False for bin_set in ref_bins_map.values()])
		freq = float(freq) / float(metadata.shape[0])
		freq_map[bid] = freq
	
	# number of bins within a reference set
	ref_n_bins = [len(ref_bins_map[acc]) for acc in metadata['accession']]
	
	# define candidate bins that are potentially imputed
	candidates = list(all_ref_bins - set(sc_bins))
	candidates.sort()
	if args.simulate:
		candidates = candidates[:100]
	candidates = np.array(candidates)
	print('Number of candidate bins: %d%s'%(len(candidates), 
	   '' if not args.simulate else ' (simulation mode is on)'))
	
	# depending on the bins in the sparse input, create the training features matrix
	training_features = np.zeros((metadata.shape[0], len(sc_bins)), dtype=np.int)
	for index, accession in enumerate(metadata['accession']):
		training_features[index] = np.array([1 if bid in ref_bins_map[accession] else 0 for bid in sc_bins])
	training_features = np.array(training_features)
	print('Shape of matrix for training features:', training_features.shape)
	
	# searching for bins with exactly the same values, for bins having the same
	# column values, only one model is trained and applied for all these bins
	bin_index_map = {}
	key_index_map = {}
	uniq_cand_labels = []
	index = 0
	for bid in candidates:
		class_labels = np.array([1 if bid in ref_bins_map[acc] else 0 for acc in metadata['accession']])
		key = ''.join(map(str, class_labels))
		if key not in key_index_map:
			uniq_cand_labels.append(class_labels)
			key_index_map[key] = index
			index += 1
		bin_index_map[bid] = key_index_map[key]
	uniq_cand_labels = np.array(uniq_cand_labels)
	
	print('\n##### Pre-Processing is done ... #####\n')
	print('Training %d models for %d candidate bins'%(uniq_cand_labels.shape[0], len(candidates)))
	print('Random Forest is used with %d trees'%(args.estimators))
	print('The task is shared by %d processors'%(size))
	
	# rank 0 shares the dimensions for other ranks used to initiate
	# the matrices containing the training features and candidates
	dimensions = (metadata.shape[0], len(sc_bins), uniq_cand_labels.shape[0])
	stdout.flush()
dimensions = comm.bcast(dimensions, root=0)
n_ref_exps, n_sc_bins, dim1_uniq_cand_labels = dimensions
if rank != 0:
	training_features = np.empty((n_ref_exps,n_sc_bins), dtype=np.int)
	uniq_cand_labels = np.empty((dim1_uniq_cand_labels, n_ref_exps), dtype=np.int)

# rank 0 pre-processed the data, created the training features, and collected the 
# candidates. This data is now shared with other ranks by broadcast
comm.Bcast(training_features, root=0)
comm.Bcast(uniq_cand_labels, root=0)

# calculate the range that defines the tasks of the rank
chunk = int(dim1_uniq_cand_labels / size) + 1
f = rank * chunk
t = (rank+1) * chunk
if t > dim1_uniq_cand_labels:
	t = dim1_uniq_cand_labels

# imp_prob_rank collects the results of each single rank
imp_prob_rank = np.zeros(chunk, dtype=np.float) - 1.0

# initialize classification algorithm with Random Forest
clf = RandomForestClassifier(n_estimators=args.estimators, random_state=42)

# for each candidate bin a classification model is trained, the more
# cores are used, the less models are trained by one rank
artificial_inst = np.array([training_features.shape[1] * [1]])
for send_index, reference_index in enumerate(range(f,t)):
	class_vector = uniq_cand_labels[reference_index]
	prob = None
	# additionally check whether the bins labels are only 0 or 1
	sum_y = sum(class_vector)
	if sum_y == 0 or sum_y == training_features.shape[0]:
		prob = 0.0 if sum_y == 0 else 1.0
	else:
		# train classification model and get the probability
		clf.fit(training_features, class_vector)
		prob = clf.predict_proba(artificial_inst)[0][1]
	imp_prob_rank[send_index] = prob

# MPI-gather of the results from each rank, collected in rank 0
if rank == 0:
	imp_prob_all = np.empty((size,chunk), dtype=np.float)
comm.Gather(imp_prob_rank, imp_prob_all, root=0)


if rank == 0:
	# extract imputed probabilities and collect reference frequencies
	index_prob = {}
	index = 0
	for r in range(size):
		for prob in imp_prob_all[r]:
			if index >= len(candidates):
				continue
			index_prob[index] = prob
			index += 1
	bin_freq_prob = []
	for bid in candidates:
		bid_index = bin_index_map[bid]
		prob = index_prob[bid_index]
		bin_freq_prob.append((bid, freq_map.get(bid,-1.0), prob))
	bin_freq_prob = sorted(bin_freq_prob, key=lambda x: x[-1], reverse=True)
	
	# add the single-cell bins to the results
	sc_bins_freq_def = [(bid, freq_map.get(bid, 0.0), -1.0) for bid in sc_bins]
	bin_freq_prob = sc_bins_freq_def + bin_freq_prob
	
	# create the SPARTAN output table and the imputed bins bed file
	impute_n = int(np.mean(ref_n_bins))
	imputed_bed = ''
	extended_bed = 'BinID\tchromosome\tstart\tend\tfrequency\timputed_probability\n'
	for n_bins_imp, (binID, frequency, imp_probability) in enumerate(bin_freq_prob):
		if n_bins_imp <= impute_n:
			imputed_bed += '%s\t%d\t%d\n'%(bin_bed_map[binID])
		
		extended_bed += str(binID) + '\t'
		extended_bed += '%s\t%d\t%d\t'%(bin_bed_map[binID])
		extended_bed += '%.5f\t%f\n'%(frequency, imp_probability)
	
	# preparing output files
	out_prefix = args.outdir
	out_prefix += '/' if out_prefix[-1] != '/' else ''
	if not os.path.exists(out_prefix):
		os.makedirs(out_prefix)
	input_file_name = args.bed.split('/')[-1].replace('.bed','')
	
	open('%s%s_imputed.bed'%(out_prefix, input_file_name), 'w').write(imputed_bed)
	open('%s%s_freq_prob.spartan'%(out_prefix, input_file_name), 'w').write(extended_bed)
	
	print('\n##### Writing output to "%s" #####\n'%(out_prefix))
	print('Reference bulk experiments have in average %d bins'%(impute_n))
	last_comment = ' (not really, because of "simulate")' if args.simulate else ''
	print('%d bins were imputed%s'%(impute_n - len(sc_bins), last_comment))
	print('Done!\n')


























