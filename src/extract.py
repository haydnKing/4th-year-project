"""Extract and clean PPRs from targets"""
import utils, targetp
from utils import pairwise
from pyHMMER import HMMER
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import copy

models = utils.loadmodels()

def extract(localization='C'):
	"""Extract all PPRs targeted to the chloroplast and clean the gaps"""
	pprs = simple_extract_all(localization)
	return pprs

def simple_extract_all(localization=None):
	files = utils.gettargetnames()
	pprs = []
	for f in files:
		pprs += simple_extract(utils.loadnuclear(f), localization)
	return pprs

def simple_extract(target, localization = None, domE=100.0):
	"""Extract all the PPRs found in target"""
	if isinstance(target, SeqRecord):
		target = [target,]

	search = HMMER.hmmsearch(hmm = models[3], targets = target, domE=domE)

	pprs = []

	for t in target:
		feats = search.getFeatures(t)
		pprs.extend(get_pprs(t, feats))	

	targetp.targetp(pprs, annotation='localization')

	if localization:
		pprs = [p for p in pprs if p.annotations['localization'] == localization]

	return clean_all(pprs)

def get_pprs(record, features):
	"""Return a list of possible PPR proteins given the target and the HMM
	matches found within it
	
		Assumes features are the same length and are made up of a single region
	"""
	#sort features into order
	features.sort(key=lambda f: int(f.location.start))

	frames = {-3: [], -2: [], -1: [], 1:[], 2:[], 3:[],}
	#split into frames
	for feature in features:
		if feature.location.strand >= 0:
			frames[int(feature.location.start) % 3 + 1].append(feature)
		else:
			frames[-((len(record) - int(feature.location.end))%3) -1].append(feature)

	#keep track of stats
	stats = {	'used': [],
			'unused': [],
			}
	
	#check for stops
	def find_stop(seq):
		"""iterate through the codons of seq"""
		for i in range(0,len(seq),3):
			if str(seq[i:i+3]).lower() in ['tga', 'tag', 'taa',]:
				return i
		return -1

	#check for stops
	def find_start(seq):
		"""iterate through the codons of seq"""
		for i in range(len(seq)-3,0,-3):
			if str(seq[i:i+3]).lower() in ['atg']:
				return i
		return -1

	def continues(old, new):
		"""Do old and new belong in the same PPR?"""
		if not old:
			return True
		if old.location.strand != new.location.strand:
			return False

		pos = (int(old.location.end), int(new.location.start))
		f = FeatureLocation(min(pos), max(pos),
				new.location.strand)

		if (pos[1] - pos[0]) > 500:
			return False
		elif find_stop(f.extract(record.seq)) >= 0:
			return False
		return True

	def chain_to_record(chain):
		"""Convert a chain of SeqFeatures to a SeqRecord"""
		#immediately discard those shorter than 2
		if len(chain) < 2:
			return None

		#extract the sequence +- margin
		pos = (int(chain[0].location.start), int(chain[-1].location.end))
		margins = (min(min(pos),1000), min(len(record) - max(pos), 1000))
		f = FeatureLocation(min(pos)-margins[0], max(pos)+margins[1], chain[0].location.strand)
		seq = f.extract(record.seq)

		#Find the start and stop codons
		start = find_start(seq[0:margins[0]])
		stop = find_stop(seq[-margins[1]:])

		#if we failed to find one
		if start < 0 or stop < 0:
			return None
		seq = seq[start:len(seq)-margins[1]+stop]

		features = []
		strand = chain[0].location.strand
		datum = (int(chain[0].location.start)-margins[0]+start if strand >= 0 else
			int(chain[-1].location.end) + margin[1] - start)

		for feature in chain:
			f = copy.deepcopy(feature)
			if(strand >= 0):
				f.location = FeatureLocation(int(f.location.start) - datum, 
					int(f.location.end) - datum, 1)
			else:
				f.location = FeatureLocation(datum - int(f.location.end),
					datum - int(f.location.start), -1)
			features.append(f)

		return SeqRecord(seq, features=features, annotations={
			'sourceid': record.id,
			'sourcestart': min(pos)-margins[0]+start,
			'sourcestrand': strand, })

	chains = []
	chain = []

	for frame,frame_features in frames.iteritems():
		#for each feature
		last_feature = None
		for feature in frame_features:
			if continues(last_feature, feature):
				chain.append(feature)
			elif chain:
				chains.append(chain)
				chain = []
			last_feature = feature
		#chains can't continue accross frames!
		if chain:
			chains.append(chain)
			chain = []

	pprs = filter(None, [chain_to_record(chain) for chain in chains])

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
		clean_gaps(ppr)
		clean_length(ppr)	

	return pprs

def extract_test():
	return simple_extract(utils.load_test())[0]

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
			if m.qualifiers.get('type', '') in stats:
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

	print "{} PPRs".format(len(pprs))	
	for i,h in enumerate(hist):
		print "{:5d}|{}".format(i+1, '*'*h)

	print ""
	print_PLS(pprs)
	
def get_ppr10():
	p = utils.load_test()
	ppr10 = simple_extract(p)
	return ppr10[0]
