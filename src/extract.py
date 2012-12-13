"""Extract and clean PPRs from targets"""
import utils, targetp
from utils import pairwise
from pyHMMER import HMMER
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

models = utils.loadmodels()

def extract(localization='C'):
	"""Extract all PPRs targeted to the chloroplast and clean the gaps"""
	pprs = simple_extract_all(localization)
	print "Found {} pprs".format(len(pprs))
	r = clean_all(pprs)
	print "Cleaning added {} extra motifs, {} gaps remain".format(*r)
	show_stats(pprs)
	return pprs

def simple_extract_all(localization=None):
	files = utils.gettargetnames()
	pprs = []
	for f in files:
		pprs += simple_extract(f, localization)
	return pprs

def simple_extract(filename, localization = None):
	sources = utils.loadnuclear(filename)
	target = utils.loadchloroplast(filename)

	search = HMMER.hmmsearch(hmm = models[3], targets = sources)
	pprs = search.getProteins(maxgap=200, mingap=-10, minlen=3, max_5_prime=1000,
			max_3_prime=1000)

	targetp.targetp(pprs, annotation='localization')

	if localization:
		t = []
		for p in pprs:
			if p.annotations['localization'] == localization:
				t.append(p)
		pprs = t

	return pprs

def get_c_terminus(pprs):
	"""Return a list of the c-termini of each protein"""
	ret = []
	if isinstance(pprs, SeqRecord):
		pprs = [pprs,]
	
	for i,ppr in enumerate(pprs):
		last = 0
		for motif in ppr.features:
			if 'PPR' in motif.type and int(motif.location.end) > last:
				last = int(motif.location.end)
		#make sure the c-terminus is in frame
		last = last - (last % 3)
		#extract and translate the c-terminus
		ret.append(SeqRecord(ppr.seq[last:].translate(), id=str(i), name="C-Terminus",
			description="C-Terminus from PPR {}".format(i)))

	if len(ret) == 1:
		return ret[0]
	return ret



def find_gaps(ppr, mingap=30, maxgap=None):
	"""Find all the gaps between PPR motifs which are gte maxgap"""
	loc = []
	feats = sorted(ppr.features, key = lambda(p): int(p.location.start))
	for a,b in pairwise(feats):
		l = abs(int(b.location.start) - int(a.location.end))

		if l >= mingap and l <= (maxgap or 'inf'):
			#We've found a gap
			g = FeatureLocation(int(a.location.end), 
				int(b.location.start), strand=1)
			g.prev = a
			g.next = b
			loc.append(g)

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
		clean_length(ppr)	
		ret = [a+b for a,b in zip(ret,r)]

	return ret

#offsets of each model compared with the actual motif start
offsets = {'PPR':-2, 'PPR_1':5, 'PPR_2':1, 'PPR_3':-1,}

def clean_length(ppr):
	"""Clean the length of the PPR motifs such that they start at the correct
	location and so that small gaps are eliminated

	adds the qualifier type as P, L or S
	"""
	
	#move each motif start position to account for the model start offset
	for f in ppr.features:
		if f.type in offsets:
			f.location = FeatureLocation(int(f.location.start)+3*offsets[f.type],
					f.location.end, f.location.strand)

	gaps = find_gaps(ppr, mingap=0, maxgap=10)
	#move the end of each motif to close small gaps between it and its successor
	for g in gaps:
		g.prev.location = FeatureLocation(g.prev.location.start,
				g.next.location.start, g.prev.strand)

	for f in ppr.features:
		if len(f) == 105: # 105 = 35 * 3
			f.qualifiers['type'] = 'P'
		elif len(f) < 105:
			f.qualifiers['type'] = 'S'
		elif len(f) > 105:
			f.qualifiers['type'] = 'L'

def print_PLS(pprs):
	stats = {'P':0,'L':0,'S':0,}

	for p in pprs:
		for m in p.features:
			if m.qualifiers['type'] in stats:
				stats[m.qualifiers['type']] += 1
	
	for k,v in stats.iteritems():
		print '{} : {}'.format(k, '#'*v)


def show_stats(pprs):
	"""Display some useful stats about the PPRs"""
	lrange = [10000, 0]
	for p in pprs:
		l = len(p.features)
		if l < lrange[0]:
			lrange[0] = len(p.features)
		if l > lrange[1]:
			lrange[1] = l
	
	hist = [0,] * (lrange[1])
	for p in pprs:
		hist[len(p.features)-1] += 1
	
	for i,h in enumerate(hist):
		print "{:5d}|{}".format(i+1, '*'*h)

	print ""
	print_PLS(pprs)
	
