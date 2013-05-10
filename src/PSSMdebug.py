"""Test the ability of a PSSM to detect a known binding site in random
sequence"""

import PSSM, extract, random

from Bio.Seq import Seq
from Bio import SeqIO
from sys import stdout
from multiprocessing import Pool, cpu_count


target = SeqIO.parse('Test_Data/PPR10_targets.fas', 'fasta').next()
ppr = extract.extract_test()

def PSSM_single_test(length=1000):

	seq = Seq(''.join(random.choice('ACGT') for x in range(length-len(target))))
	pos = random.randrange(0, len(seq))
	seq = seq[0:pos] + str(target.seq) + seq[pos:]

	aln = PSSM.search(ppr,seq,gaps=2)
	rank = -1
	for i,a in enumerate(aln):
		if a.pos == pos+1:
			rank = i
			break
	
	return rank

def PSSM_sample(length, samples=20, pool=None):
	if not pool:
		pool = Pool(processes=cpu_count())

	it = pool.imap(PSSM_single_test, (length,) * samples)
	r = []
	for out in it:
		r.append(out)
		stdout.write("\r{:.1f}%    ".format(100.0 * float(len(r)) / float(samples)))
		stdout.flush()
	stdout.write("\r         \r")
	stdout.flush()

	return r

def write_stats(fname, data):
	keys = data.keys()
	l = len(data[keys[0]])
	for key in keys:
		if len(data[key]) != l:
			raise ValueError("Values must be the same length")

	f = open(fname, 'w')
	#write header line
	f.write(("{:>10s} "*len(keys) + '\n').format(*keys))

	for i in range(l):
		for key in keys:
			f.write("{:10.5f} ".format(data[key][i]))
		f.write('\n')

	f.close()

def test(samples=50):
	import matplotlib.pyplot as plt
	
	lengths = range(1000,100000,10000)

	pool = Pool(processes=cpu_count())
	correct = []
	top10 = []
	missed = []

	for l in lengths:
		print "l = {}".format(l)
		ret = PSSM_sample(l,samples,pool)
		c = 0
		t = 0
		m = 0
		for r in ret:
			if r < 0:
				m += 1
			elif r < 10:
				t += 1
			if r == 0:
				c += 1
		correct.append(float(c) / float(samples))
		top10.append(float(t) / float(samples))
		missed.append(float(m) / float(samples))

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(	lengths, correct, 'g', 
						lengths, top10, 'b', 
						lengths, missed, 'r')
	ax.legend(['Correct', 'Top 10', 'Not Found'], 'upper right')
	ax.set_ylim([-0.1,1.1])
	ax.set_ylabel('Probability')
	ax.set_xlabel('Target Length (bp)')
	plt.savefig('output/PSSMsearch.png')
	plt.savefig('output/PSSMsearch.svg')

	write_stats('output/PSSMsearch.csv', {	
		"length": lengths, 
		"correct": correct, 
		"top10": top10, 
		"missed": missed,
	})
	
if __name__ == '__main__':
	test()
