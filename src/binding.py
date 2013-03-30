"""Attempt to match PPRs based on their binding preferences"""

import extract, utils, PSSM
import ppr as PPR

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

def show_distance(num=0):
	genomes = PPR.load_genomes()

	ppr = genomes[0].pprs[num]
	closest = minima(ppr, genomes)
	left = range(len(closest[0]))

	plt.bar(left, closest[1], width=1.0)
	plt.xticks([l+0.5 for l in left], closest[0], rotation=90)
	plt.show()

def ppr_distance():
	genomes = PPR.load_genomes()

	data = []

	for i,g in enumerate(genomes, 1):
		print "Genome {}/{}: {} ({})".format(i, len(genomes),
				g.name,reduce_name(g.name))
		data.append([reduce_name(g.name),] + get_average(g,genomes))


	utils.write_data(['',] + [reduce_name(g.name) for g in genomes],
			data,
			'output/Average_PPR_distance.dat')

def get_groups(threshold=0.3):

	stdout.write("\rLoading...")
	stdout.flush()

	genomes = PPR.load_genomes()

	groups = []

	total = sum([sum([1 for p in g.pprs]) for g in genomes])
	done = 0
	fmt = "\rGrouping ({{:{}d}}/{{}}) | {{:2.1f}}%, {{:.1f}} per group     ".format(
			len(str(total)))
	for genome in genomes:
		for ppr in genome.pprs:
			#find the best group
			possible = []
			#groups are possible if the distance to each member is lte threshold
			for group in groups:
				if max([ppr.distance(p2, 'MSE')[0] for p2 in group]) <= threshold:
					possible.append(group)
			#if there are no possible groups, start a new one
			if not possible:
				groups.append([ppr,])
			#otherwise, append to the group with the lowest average distance
			else:
				avg = [sum([ppr.distance(p2, 'MSE')[0] for p2 in group])/float(len(group))
									for group in possible]
				possible[avg.index(min(avg))].append(ppr)
			#Write progress info
			done += 1
			stdout.write(fmt.format(done,total, 100.0 * done / float(total),
				done / float(len(groups))))
			stdout.flush()

	stdout.write("\r" + (" "*60) + "\r")
	stdout.flush()
	return groups

def group_stats(groups):
	l = [len(g) for g in groups]
	print "{} groups, {} <= length <= {}, average = {}".format(
			len(groups), min(l), max(l), sum(l)/float(len(l)))

if __name__ == '__main__':
	ppr_distance()


