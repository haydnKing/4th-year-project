"""Attempt to match PPRs based on their binding preferences"""

import extract, utils, ppr, PSSM

import matplotlib.pyplot as plt

from sys import stdout

def minima(ppr, genomes):
	names = []
	closest = []
	for g in genomes:
		c = g.closest(ppr)
		names.append(g.name.replace('_', '\n'))
		closest.append(ppr.distance(c)[0])
	return (names, closest)

def get_average(genome, genomes):
	total = [0,]*len(genomes)

	fmt = "\r({{:>{}d}}/{{}}) {{:3.0f}}%".format(
			len(str(len(genome.pprs))))

	for i,p in enumerate(genome.pprs, 1):
		stdout.write(fmt.format(i,len(genome.pprs),
			100.0*i/float(len(genome.pprs))))
		stdout.flush()
		names, closest = minima(p, genomes)
		total = [c+t for c,t in zip(closest,total)]

	stdout.write("\r" + (" "*30) + "\r")
	stdout.flush()
	return [t/len(genome.pprs) for t in total]

def reduce_name(name):
	return "".join([s[0].upper() for s in name.split('_')])

def ppr_distance():
	genomes = ppr.load_genomes()

	data = []

	for i,g in enumerate(genomes, 1):
		print "Genome {}/{}: {} ({})".format(i, len(genomes),
				g.name,reduce_name(g.name))
		data.append([reduce_name(g.name),] + get_average(g,genomes))


	utils.write_data(['',] + [reduce_name(g.name) for g in genomes],
			data,
			'output/Average_PPR_distance.dat')


if __name__ == '__main__':
	ppr_distance()


