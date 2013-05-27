"""Design a PSSM matrix and search for it in the target sequence"""

from math import log, exp, sqrt, pi
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import matplotlib.pyplot as plt
from sys import stdout

import binding_rules

class Alignment:

	def __init__(self, pos, odds, gaps):
		self.pos = pos
		self.odds = odds
		self.gaps = gaps

	def toSeqFeature(self, pssm, reflect=None):
		l = len(pssm) + len(self.gaps)
		if not reflect:
			return SeqFeature(FeatureLocation(self.pos, self.pos+l, strand=1), 
					type="PSSM", qualifiers={"odds":self.odds})
		else:
			return SeqFeature(FeatureLocation(reflect-self.pos-l, reflect-self.pos,
				strand=-1),
				type = "PSSM", qualifiers={"odds":self.odds})

	def __lt__(self, other):
		return self.odds < other.odds

	def __le__(self, other):
		return self.odds <= other.odds

	def __gt__(self, other):
		return self.odds > other.odds

	def __ge__(self, other):
		return self.odds >= other.odds

def search(ppr, target, reverse=True, gaps=1, show_stats=False):
	itarget = get_iseq(target)
	bg = get_background(itarget)

	if isinstance(ppr, SeqRecord):
		pssm = binding_rules.build_PSSM(ppr, coding='yagi')#, background=bg)
	else:
		#assume we were passed a pssm
		pssm = ppr

	alignments = PSSM_gapped_search(pssm, itarget, gaps)
	ralignments = PSSM_gapped_search(pssm, ireverse_complement(itarget), gaps)

	sf  = [x.toSeqFeature(pssm) for x in alignments]
	sf += [x.toSeqFeature(pssm,reflect=len(itarget)) for x in ralignments]
	sf.sort(key=lambda f: -f.qualifiers['odds'])

	if show_stats:
		display_results(pssm, sf, len(itarget))

	return sf

def build(sr, coding='barkan'):
	return binding_rules.build_PSSM(sr, coding)

def as_string(pssm):
	return "\n".join([("{:3.1f} "*4).format(*i) for i in pssm])

def align(a, b):
	dist, align = distance(a,b,'MSE')
	print "Dist,Align = {},{}".format(dist,align)
	if len(a) > len(b):
		l_offset = 0
		r_offset = len(a)-len(b)
	else:
		l_offset = len(b)-len(a)
		r_offset = 0

	for i in range(max(len(a),len(b))):
		l = i-l_offset
		r = i-r_offset
		if l in range(len(a)):
			line = ("{:3.1f} "*4).format(*a[l])
		else:
			line = "    "*4
		line += "| "
		if r in range(len(b)):
			line += ("{:3.1f} "*4).format(*b[r])
		print line

	print "Distance: {}".format(dist)

def display_alignments(pssm, itarget, aln):
	for i,a in enumerate(aln):
		print "Alignment {}".format(i)
		print str_alignment(pssm, itarget, a)

def str_alignment(pssm, iseq, aln):
	pssm = PSSM_insert_gaps(pssm,aln.gaps)	
	iseq = iseq[aln.pos:aln.pos+len(pssm)]
	odds = PSSM_match_odds(pssm,iseq)
	return ("Position {}: odds = {} gaps = {}\n\tconsensus: {}\n\t           {}"+
		"\n\ttarget   : {}").format(
			aln.pos,aln.odds,aln.gaps,
			PSSM_consensus_seq(pssm),
			"".join(['+' if i > 0 else '-' for i in odds]), 
			get_seq(iseq))
	#odds.sort(key=lambda o: -o[1])

def display_results(pssm, features, target_length):
	"""Show histograms of the output"""
	def linspace(low, high, num):
		step = float(high - low) / float(num)
		space = [0]*num
		for i in range(num):
			space[i] = low + step*float(i)
		return space
	
	r = [f.qualifiers['odds'] for f in features]
	(mean, var) = profile(pssm, runs=3, length=target_length)
	x = linspace(min(r), max(r), 100)
	y = [(1 / (sqrt(2 * pi * var))) * 
				exp(-0.5 * (_x - mean)*(_x-mean) / var) for _x in x]

	plt.clf()
	plt.hold(True)
	plt.plot(x,y)
	plt.hist(r, 20, normed=True)
	plt.show()

def KL_divergence(P, Q):
	"""Calculate the minimum KL divergence of to PSSMs"""
	if len(P) != len(Q):
		return None
	return sum([
		sum([(p_ij - q_ij)*exp(p_ij) for p_ij,q_ij in zip(p_i,q_i)]) 
			for p_i,q_i in zip(P,Q)
	])

def MSE(P,Q):
	"""Return the MSE of P and Q"""
	if len(P) != len(Q):
		return None
	return sum([
		sum([(p_ij - q_ij)**2 for p_ij,q_ij in zip(p_i,q_i)]) 
			for p_i,q_i in zip(P,Q)
	]) / (4.0 * len(P))

def distance(P,Q,method='KL'):
	methods = {'KL': KL_divergence,
							'MSE': MSE,}
	div = methods[method]

	if len(P) == len(Q):
		return (div(P,Q), 0)
	ret = []
	pad = abs(len(P) - len(Q))
	eq = [(log(0.25),)*4,]
	if len(P) > len(Q):
		for i in range(pad+1):
			ret.append(div(P, eq*i + Q + eq*(pad-i)))
	if len(P) < len(Q):
		for i in range(pad):
			ret.append(div(eq*i + P + eq*(pad-i), Q))
	return (min(ret), ret.index(min(ret)),)

def PSSM_insert_gaps(pssm, gaps):
	rpssm = [i for i in pssm]
	for g in gaps:
		rpssm = rpssm[0:g] + [gap,] + rpssm[g:]
	return rpssm

def PSSM_consensus_iseq(PSSM):
	return [i.index(max(i)) for i in PSSM]

def PSSM_consensus_seq(PSSM):
	"""Return the consensus -
		UPPER: largest odds
		lower: closely tied
		dash-: equal odds
	"""
	ret = ""
	for s in PSSM:
		b = int2base[s.index(max(s))]
		s = list(reversed(sorted(s)))
		#if s is a gap
		if s == [0,0,0,0,]:
			ret += '.'
		#if s is equal probability
		elif (s[0] - s[3]) <= 0.001:
			ret += '-'
		#if s shows a strong preference
		elif (s[0] - s[1]) > 0.5:
			ret += b.upper()
		#else s shows a weak preference
		else:
			ret += b.lower()
	return ret

def PSSM_match_odds(PSSM, iseq):
	return [p[b] for p,b in zip(PSSM,iseq)]

def PSSM_gapped_search(PSSM, itarget, max_gaps=5):
	#do an exhaustive search, and store the maxima
	maxima = non_max_suppression(PSSM_exhaustive_search(PSSM,itarget))
	l = len(maxima)

	alignments = []

	#for each maximum
	for i,aln in enumerate(maxima):
		#extract the envelope
		f = max(0, aln.pos - max_gaps)
		t = min(len(itarget), aln.pos + len(PSSM) + max_gaps)
		iseq = itarget[f:t+1]

		#localise within the envelope
		aln = PSSM_localise(PSSM, iseq, max_gaps)
		aln.pos += f
		alignments.append(aln)

	return list(reversed(sorted(alignments)))

def PSSM_localise(PSSM, iseq, max_gaps):
	"""Find the best alignment within the sequence"""
	if len(PSSM) > len(iseq) or max_gaps < 0:
		return None
	#Find my best alignment
	aln = PSSM_find_max(PSSM,iseq)

	if max_gaps == 0 or not aln:
		return aln

	for g in range(1,len(PSSM)):
		aln2 = PSSM_localise(PSSM[0:g] + [gap,] + PSSM[g:], iseq, max_gaps-1)
		if not aln2:
			continue
		if aln2 > aln:
			aln2.gaps = [g,] + aln2.gaps
			aln = aln2

	return aln

def PSSM_find_max(PSSM, itarget):
	r = PSSM_exhaustive_search(PSSM, itarget)
	if r:
		s = max(r)
		return Alignment(r.index(s), s, [])
	return None

def PSSM_exhaustive_search(PSSM, itarget):
	ret = []
	for i in range(0, len(itarget) - len(PSSM)):
		ret.append(sum([p[b] for p,b in zip(PSSM, itarget[i:i+len(PSSM)])]))
	return ret

def non_max_suppression(odds):
	if len(odds) == 1:
		return [Alignment(0,odds[0],[]),]
	ret = []
	def is_max(i):
		return odds[i] > odds[i-1] and odds[i] > odds[(i+1)%len(odds)]

	return [Alignment(i,o,[]) for i,o in enumerate(odds) if is_max(i)]

def get_background(itarget):
	"""calculate background probabilities"""
	#add a prior for calculating the background with minimal data
	background = [10,10,10,10]
	for i in itarget:
		background[i] += 1
	return tuple(float(i) / sum(background) for i in background)

def get_iseq(seq):
	if isinstance(seq, SeqRecord):
		return [base2int[i] for i in str(seq.seq).upper()]
	return [base2int[i] for i in str(seq).upper()]

def get_seq(iseq):
	return "".join((int2base[i] for i in iseq))

def profile(pssm, runs=3, gaps=1, background=[1,1,1,1], length=None):
	out = []
	for i in range(runs):
		stdout.write('\rProfile: {} / {}   '.format(i+1, runs))
		stdout.flush()
		out += _random(pssm,gaps,background,length)
	stdout.write('\r                          \r')
	stdout.flush()
	mean = sum(out) / float(len(out))
	var  = sum([ (o-mean)*(o-mean) for o in out]) / float(len(out))

	x = range(int(min(out)), int(max(out))+1)
	y = [(1 / (sqrt(2*pi*var))) * 
				exp(-0.5 * (_x-mean)*(_x-mean) / var) for _x in x]
	plt.clf()
	plt.hold(True)
	plt.hist(out, bins=20, normed=True)
	plt.plot(x,y)
	plt.show()
	return (mean, var)


import random
def _random(pssm, gaps, background, l=None):
	if not l:
		l = len(pssm) + gaps
	b = [x/float(sum(background)) for x in background]

	#generate the sequence according to the background statistics
	iseq = [0] * l
	for i in range(l):
		r = random.random()
		for j in range(4):
			if r < sum(b[0:j+1]):
				break
		iseq[i] = j

	alignments = PSSM_gapped_search(pssm, iseq, gaps)
	alignments += PSSM_gapped_search(pssm, ireverse_complement(iseq), gaps)

	return [x.odds for x in alignments]

def ireverse_complement(iseq):
	return [icomplement[b] for b in reversed(iseq)]

base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3,}
int2base = ('A','C','G','T')
icomplement = (3, 2, 1, 0)

equal = (1, 1, 1, 1)
log_equal = tuple([log(0.25),]*4)
gap = (0,0,0,0,)

