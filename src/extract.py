"""Extract and clean PPRs from targets"""
import utils, targetp
from pyHMMER import HMMER
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

models = utils.loadmodels()

def extract():
	"""Extract all PPRs targeted to the chloroplast and clean the gaps"""
	pprs = simple_extract_all()
	print "Before Clean"
	show_stats(pprs)
	r = clean_all(pprs)
	print "Cleaning added {} extra motifs, {} gaps remain".format(*r)
	show_stats(pprs)
	return pprs

def simple_extract_all(localization='C'):
	files = utils.gettargetnames()
	pprs = []
	for f in files:
		pprs += simple_extract(f, localization)
	return pprs

def simple_extract(filename, localization = None):
	sources = utils.loadnuclear(filename)
	target = utils.loadchloroplast(filename)

	search = HMMER.hmmsearch(hmm = models[0], targets = sources)
	pprs = search.getProteins(maxgap=100, mingap=0, minlen=4, max_5_prime=500,
			max_3_prime=500)

	targetp.targetp(pprs, annotation='localization')

	if localization:
		t = []
		for p in pprs:
			if p.annotations['localization'] == localization:
				t.append(p)
		pprs = t

	return pprs

import itertools
def pairwise(iterable):
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = itertools.tee(iterable)
	next(b, None)
	return itertools.izip(a, b)

def find_gaps(ppr, maxgap=30):
	"""Find all the gaps between PPR motifs which are gte maxgap"""
	loc = []
	feats = sorted(ppr.features, key = lambda(p): int(p.location.start))
	for a,b in pairwise(feats):
		if (int(b.location.start) - int(a.location.end)) >= maxgap:
			#We've found a gap
			loc.append(FeatureLocation(int(a.location.end), 
				int(b.location.start), strand=1))

	return loc

def clean_gaps(ppr):
	"""find and remove as many gaps as possible"""
	ret = [0, 0]
	repeat = True
	while repeat:
		repeat = False
		gaps = find_gaps(ppr)
		#search for matches to the other models in each gap
		for g in gaps:
			r = SeqRecord(g.extract(ppr.seq))
			search = HMMER.hmmsearch(hmm=models, targets=r)
			for m in search.matches:
				if m.frame >= 0:
					if not m.withinTarget():
						continue
					repeat = True
					ret[0] = ret[0] + 1
					#get the feature which refers to the match, offset to the ppr
					#and add it to the ppr
					ppr.features.append(m.asSeqFeature(mode='hmm', offset=-g.start))
					break
	ret[1] = len(gaps)
	ppr.features = sorted(ppr.features, key = lambda(p): int(p.location.start))
	return ret

def clean_all(pprs):
	ret = [0,0]
	for ppr in pprs:
		r = clean_gaps(ppr)
		ret = [a+b for a,b in zip(ret,r)]

	return ret

def show_stats(pprs):
	"""Display some useful stats about the PPRs"""
	print "{} PPRs".format(len(pprs))
	lrange = [10000, 0]
	for p in pprs:
		l = len(p.features)
		if l < lrange[0]:
			lrange[0] = len(p.features)
		if l > lrange[1]:
			lrange[1] = l
	
	print "Shortest: {}, Longest: {}".format(*lrange)
	hist = [0,] * (lrange[1])
	for p in pprs:
		hist[len(p.features)-1] += 1
	
	for i,h in enumerate(hist):
		print "{:5d}|{}".format(i+1, '*'*h)
	
