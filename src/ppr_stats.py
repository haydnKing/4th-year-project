import utils, os, os.path
from Bio import SeqIO


def get_genome_names():
	genomes = os.listdir('output/PPRs/')
	return sorted([os.path.splitext(f)[0] for f in genomes])

def load_pprs(genome):
	path = "output/PPRs/{}.gb".format(genome)
	return list(SeqIO.parse(path, "genbank"))

def length_hist(pprs):
	if not pprs:
		return [0,]
	lengths = [len(ppr.features) for ppr in pprs]
	hist = [0,] * (1+max(lengths))
	for l in lengths:
		hist[l] += 1
	return hist

def fmt_hist(hists):
	max_length = max([len(h) for h in hists]) - 1
	data = []
	for l in range(2, max_length+1):
		row = [l,]
		for h in hists:
			try:
				row.append(h[l])
			except IndexError:
				row.append(0)
		data.append(tuple(row))
	return data

def get_stats(pprs):
	return (len(pprs),)

def get_stat_header():
	return ("number",)

def main():
	genomes = get_genome_names()
	hist = []
	data = []
	for g in genomes:
		pprs = load_pprs(g)
		data.append((g,) + get_stats(pprs))
		hist.append(length_hist(pprs))

	hist = fmt_hist(hist)

	utils.write_data(('genome',)+get_stat_header(), data, "output/ppr_data.dat")
	utils.write_data(['length',]+genomes, hist, "output/ppr_lengths.dat")

if __name__ == "__main__":
	main()
