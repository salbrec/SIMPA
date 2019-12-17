
def get_peaks(file_path, allowed_chromosomes, enrich_index=8):
	allowed_chromosomes_set = set(allowed_chromosomes)
	peaks = {}
	for chrom in allowed_chromosomes:
		peaks[chrom] = []
	with open(file_path, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) < 3:
				return None
			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			fdr = -1
			if enrich_index != -1:
				fdr = float(line[enrich_index])
			if chrom in allowed_chromosomes_set:
				peaks[chrom].append((chrom, start, end, fdr))
	return peaks

def bin_it(peaks, allowed_chromosomes, chrom_sizes, bin_size):
	# create each bin and check whether a peak overlaps it
	bins_set = set()
	position_map = {}
	bin_fdr = []
	bin_ID = 0
	for chrom in allowed_chromosomes:
		chrom_size = chrom_sizes[chrom]
		bin_end = bin_size
		while bin_end < chrom_size:
			bin_start = bin_end - bin_size
			for peak in peaks[chrom]:
				if peak[2] > bin_start:
					if peak[1] < bin_end:
						bins_set.add(bin_ID)
						bin_fdr.append((bin_ID, peak[3]))
					else:				# early exit: since peaks are sorted the loop can be broken
						break			#             when the start of the bin is "behind" the end of the current peak
			output_end = bin_end if bin_end < chrom_size else chrom_size
			position_map[bin_ID] = (chrom, bin_start, output_end)
			# update the running values
			bin_end += bin_size
			bin_ID += 1
	return bins_set, bin_fdr, bin_ID, position_map

