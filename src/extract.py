#!/usr/bin/python

"""Extract and clean PPRs from targets"""
import utils, targetp, ppr, classify
from utils import pairwise
from pyHMMER import HMMER
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from Bio import SeqIO

import os.path
import copy
import sys

models = utils.loadmodels()

def extract(localization=None, files=None, verbose=False):
	"""Extract all PPRs targeted to the chloroplast and clean the gaps"""
	pprs = simple_extract_all(localization, files, verbose)
	return pprs

def simple_extract_all(localization=None, files=None, verbose=False):
	if not files:
		files = utils.gettargetnames()
	pprs = []
	for f in files:
		if verbose:
			print "Processing {}".format(f)
		tpprs = []
		for i,c in enumerate(utils.loadnuclear(f)):
			if verbose:
				print "\t{} ({:,d} bp)".format(i,len(c))
			tpprs += simple_extract(c, localization)
		genome = os.path.splitext(os.path.basename(f))[0]
		tpprs.sort(key=lambda p: -len(p.features))
		for i,p in enumerate(tpprs):
			p.id = "PPR{}|{}".format(i,genome)
		pprs += tpprs
		if verbose:
			print "\tFound {}, {} total".format(len(tpprs), len(pprs))
	pprs.sort(key=lambda p: -len(p.features))
	return pprs

def annotate_pprs(target, remove = False, localization = None):
	"""Annotate the PPRs found in the target"""
	print "Searching..."
	#find all easy-to-locate PPR motifs
	search = HMMER.hmmsearch(hmm = models[3], targets = target)

	#get features for each motif
	motifs = search.getFeatures(target)

	print "Got {} motifs, grouping...".format(len(motifs))
	#group features by frame and locatiion
	groups = group_motifs(motifs, max_gap=1000)

	print "Got {} groups, extracting envelopes...".format(len(groups))
	#extract the sequence envelope around each group
	envelopes = [get_envelope(group, target, margin=500) for group in groups]

	print "Got {} envelopes, locating PPRs...".format(len(envelopes))
	#locate the PPR within each group
	pprs = [locate_ppr(envelope) for envelope in envelopes]

	feats = [SeqFeature(FeatureLocation(
			p.annotations['src_from'], 
			p.annotations['src_to'],
			p.annotations['src_strand']),
			"PPR_protein",
			id="PPR_{}".format(i)) for i,p in enumerate(pprs)]
	if remove:
		target.features = []
	target.features += feats
	return target

def simple_extract(target, localization = None):
	"""Extract all the PPRs found in target"""
	if not isinstance(target, SeqRecord):
		raise TypeError("simple_extract requires a Bio.SeqRecord, not {}".format(
			type(target)))

	print "Searching..."
	#find all easy-to-locate PPR motifs
	search = HMMER.hmmsearch(hmm = models[3], targets = target)

	#get features for each motif
	motifs = search.getFeatures(target)

	print "Got {} motifs, grouping...".format(len(motifs))
	#group features by frame and locatiion
	groups = group_motifs(motifs, max_gap=1000)

	print "Got {} groups, extracting envelopes...".format(len(groups))
	#extract the sequence envelope around each group
	envelopes = [get_envelope(group, target, margin=1000) for group in groups]

	print "Got {} envelopes, locating PPRs...".format(len(envelopes))
	#locate the PPR within each group
	pprs = [locate_ppr(envelope) for envelope in envelopes]
	pprs = [p for p in pprs if p != None]
	pprs = [add_source(p, target) for p in pprs]

	print "Got {} PPRs, cleaning...".format(len(groups))
	#clean the gaps between features
	pprs = [clean_gaps(ppr) for ppr in pprs]

	#annotate the tail region and classify each PPR
	classify.classify(pprs)

	#predict each PPR's target
	targetp.targetp(pprs, annotation='localization')

	#filter the desired location
	if localization:
		pprs = [p for p in pprs if p.annotations['localization'] == localization]

	#DEBUG
	feats = [SeqFeature(FeatureLocation(
			p.annotations['src_from'], 
			p.annotations['src_to'],
			p.annotations['src_strand']),
			"envelope",
			qualifiers = {
				'name':"env_{}".format(i),
			}) 
				for i,p in enumerate(envelopes)]
	feats += [SeqFeature(FeatureLocation(
			p.annotations['src_from'], 
			p.annotations['src_to'],
			p.annotations['src_strand']),
			"PPR",
			qualifiers = {
				'name':"ppr_{}".format(i),
			}) 
				for i,p in enumerate(pprs)]
	feats += motifs
	target.features = feats
	print len(target.features)
	SeqIO.write(target, "simple_extract.gb", 'genbank')

	#return a list of nicely presented PPRs
	return pprs

def group_motifs(motifs, max_gap):
	"""Given a list of PPR motifs, group them by frame and then by location such
	the the maximum gap within a group is max_gap"""
	def distance(a,b):
		"""return the minimum distance between to seqfeatures"""
		return min(	abs(a.location.start - b.location.end),
								abs(b.location.start - a.location.end))

	by_frame = {1: [], 2: [], 3: [], -1: [], -2: [], -3: [],}
	for motif in motifs:
		by_frame[motif.qualifiers['frame']].append(motif)

	groups = []
	#for each frame
	for frame in by_frame.itervalues():
		#sort by start location
		frame.sort(key=lambda m: m.location.start)
		#while there are still motifs to be grouped
		while frame:
			#begin a new group
			current_group = [frame[0]]
			#find each motif within range
			for motif in frame[1:]:
				#we put these in ascending order, so the nearest is the last
				if distance(motif,current_group[-1]) < max_gap:
					current_group.append(motif)
			#remove each item in the current group from the frame
			for m in current_group:
				frame.remove(m)
			#add current group to the list of groups
			groups.append(current_group)

	return groups

def get_envelope(group, target, margin):
	"""Get the sequence surrounding the group"""
	#margin must be divisible by 3
	margin = margin - margin % 3

	#make sure we're still ordered
	group.sort(key=lambda k: k.location.start)
	start = max(0, group[0].location.start - margin)
	end   = min(len(target), group[-1].location.end + margin)
	frame = group[0].qualifiers['frame']

	seq = target.seq[start:end]
	if frame < 0:
		seq = seq.reverse_complement()
	seq.alphabet=generic_dna

	return SeqRecord(seq, id='env', name='PPR_env', description='PPR_envelope',
			annotations = {'src_from': start, 'src_to': end, 
				'src_strand': 1 if (frame > 0) else -1})

def locate_ppr(envelope):
	"""Find and annotate the protein within"""
	#find all the PPR motifs
	search = HMMER.hmmsearch(hmm = models[3], targets = envelope)
	motifs = search.getFeatures(envelope)

	#We're only interested in what's in frame 1 - that's what the envelope is
	#centered on
	frame = 1

	#check if there were any other frames, and drop them
	l = len(motifs)
	motifs = [m for m in motifs if m.qualifiers['frame']==frame]
	if len(motifs) < l:
		print ("WARNING: dropping some PPR motifs after finding multiple frames " +
				 "in the same envelope.")

	#A ppr must contain 2 or more PPR motifs
	if len(motifs) < 2:
		return None
	
	#order the motifs
	motifs.sort(key=lambda m: m.location.start)
	#find start codon
	start = motifs[0].location.start
	while start > 0 and str(envelope.seq[start:start+3]).lower() != "atg":
		start -= 3
	if start < 0:
		start = 0

	#find stop codon
	stop = motifs[-1].location.end
	while stop < len(envelope) and (
		str(envelope.seq[stop:stop+3]).lower() not in ["tag", "tga", "taa"]):
		stop += 3
	if stop > len(envelope):
		stop = len(envelope)

	#move the motifs
	for m in motifs:
		m.location = FeatureLocation(m.location.start-start, m.location.end-start)
	

	#get absolute start and end
	if envelope.annotations['src_strand'] > 0:
		src_from = envelope.annotations['src_from'] + start
		src_to   = envelope.annotations['src_from'] + stop
	else:
		src_from = envelope.annotations['src_to'] - stop
		src_to   = envelope.annotations['src_to'] - start

	#return a record
	return SeqRecord(envelope.seq[start:stop],
			features = motifs,
			annotations = {
				"src_from"	: src_from,
				"src_to"		: src_to,
				"src_strand": envelope.annotations['src_strand'],
			})

def add_source(ppr, source):
	"""Add source metadata to the ppr"""
	ppr.annotations.update({
		'src_id': source.id,
		'src_name': source.name,
		'src_desc': source.description,
		})
	ppr.id = "PPR"
	ppr.name = "PPR Protein"
	ppr.description = "PPR Protein"
	return ppr

#offsets of each model compared with the actual motif start
offsets = {'PPR':-2, 'PPR_1':5, 'PPR_2':1, 'PPR_3':-1,}

def clean_gaps(ppr):
	"""Clean the length of the PPR motifs such that they start at the correct
	location and so that small gaps are eliminated

	adds the qualifier type as P, L or S
	"""

	#move each motif start position to account for the model start offset
	for f in ppr.features:
		f.location = FeatureLocation(int(f.location.start)+3*offsets['PPR_3'],
			f.location.end, f.location.strand)

	gaps = find_gaps(ppr, mingap=0, maxgap=10)
	#move the end of each motif to close small gaps between it and its successor
	for g in gaps:
		g.prev.location = FeatureLocation(g.prev.location.start,
			g.next.location.start,
			g.prev.strand)

	for f in ppr.features:
		if len(f) == 105: # 105 = 35 *3
			f.qualifiers['type'] = 'P'
		elif len(f) < 105:
			f.qualifiers['type'] = 'S'
		elif len(f) > 105:
			f.qualifiers['type'] = 'L'
	
	return ppr

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

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print "Usage: {} file1 [file2...]".format(sys.argv[0])
		sys.exit(-1)
	
	files = sys.argv[1:]
	for f in files:
		pprs = extract(files=[f,], verbose = True)

		head,tail = os.path.split(f)
		genome = os.path.splitext(tail)[0]
		ppr.write_pprs(pprs, genome)

		print "###########################################################"
		print "Genome \'{}\' summary:".format(genome)
		print "-----------------------------------------------------------"
		show_stats(pprs)
		print "###########################################################"
		

