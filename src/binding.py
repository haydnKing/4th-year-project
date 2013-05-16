"""Attempt to match PPRs based on their binding preferences"""

import extract, utils, PSSM
import ppr as PPR

from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import matplotlib.pyplot as plt

import os.path
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
	"""Attempt to group similar PPRs"""

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

def annotate_all(genome_name='Arabidopsis_thaliana',plastid_name='NC_000932'):
	pprs = [p for p in PPR.load_records(genome_name) if (len(p.features) >= 15 
					and p.annotations.get('localization','').upper() == 'C')]


	plastid = utils.load_plastid(plastid_name)

	plastid.features = [f for f in plastid.features if f.type.lower() == 'gene']

	annotate_binding_domains(pprs, plastid)

	ofile = os.path.join(utils.OutDir, "Binding/{}--{}.gb".format(genome_name,
		plastid_name))

	SeqIO.write(plastid, ofile, 'gb')

def annotate_binding_domains(pprs, plastid):
	for i,ppr in enumerate(pprs):
		stdout.write("\r{}/{}   ".format(i,len(pprs)))
		stdout.flush()
		feats = get_domains(ppr, plastid, 10.0, 1)
		for f in feats:
			f.type = "PPR_{}".format(i)
		plastid.features += feats
	stdout.write("\r                          \r")
	stdout.flush()
		
def get_domains(ppr, plastid, percentile=10.0, gaps=1):
	"""Get a list of putative binding domains in the plastid"""
	a = PSSM.search(ppr,plastid,gaps=gaps)
	b = PSSM.search(ppr,plastid.reverse_complement(),gaps=gaps)
	for f in b:
		f.location = FeatureLocation(len(plastid)-f.location.end,
																 len(plastid)-f.location.start, strand=-1)
	
	feats = a + b
	feats.sort(key=lambda c: -c.qualifiers['odds'])

	top = feats[0].qualifiers['odds']
	bottom = feats[-1].qualifiers['odds']
	threshold = top - (percentile/100.0) * (top-bottom)
	
	for i in range(len(feats)):
		if feats[i].qualifiers['odds'] < threshold:
			return feats[0:i-1]

	return feats

def write_domains(doms, fname):
	f = open(fname, 'w')
	f.write("name, start, end, strand\n")

	for d in doms:
		f.write("{}, {}, {}, {}\n".format(d.type, d.location.start,
			d.location.end, d.location.strand))
	
	f.close()

if __name__ == '__main__':
	ppr_distance()


