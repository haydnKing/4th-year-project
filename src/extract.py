#!/usr/bin/python

"""Extract and clean PPRs from targets"""
import utils, targetp, ppr, classify
from utils import pairwise
from pyHMMER import HMMER
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from Bio import SeqIO, Entrez, Alphabet
Entrez.email = 'hjk38@cam.ac.uk'

import os.path
import copy
import sys

models = utils.loadmodels()

def from_entrez(gid, db='nuccore'):
	hnd = Entrez.efetch(db=db, id=gid, rettype='gb', retmode='text')
	seq = SeqIO.read(hnd, 'genbank', alphabet=Alphabet.generic_dna)
	return simple_extract(seq)

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
	pprs = simple_extract(target, localization)

	if remove:
		target.features = []

	target.features += [SeqFeature(FeatureLocation(
			p.annotations['src_from'], 
			p.annotations['src_to'],
			p.annotations['src_strand']),
			"PPR_protein",
			qualifiers = {
				'name':"ppr_{}".format(i),
			}) 
				for i,p in enumerate(pprs)]
	for p in pprs:
		target.features += [SeqFeature(FeatureLocation(
			p.annotations['src_from'] + m.location.start
				if p.annotations['src_strand'] > 0 else
			p.annotations['src_to'] - m.location.end, 
			p.annotations['src_from'] + m.location.end
				if p.annotations['src_strand'] > 0 else
			p.annotations['src_to'] - m.location.start,
			p.annotations['src_strand']),
			"PPR_Motif") 		
					for m in p.features]

	return target

def simple_extract(target, localization = None, verbose=False):
	"""Extract all the PPRs found in target"""
	if not isinstance(target, SeqRecord):
		raise TypeError("simple_extract requires a Bio.SeqRecord, not {}".format(
			type(target)))

	if verbose:
		print "Searching..."
	#find all easy-to-locate PPR motifs
	search = HMMER.hmmsearch(hmm = models[3], targets = target)

	#get features for each motif
	motifs = search.getFeatures(target)

	if verbose:
		print "Got {} motifs, grouping...".format(len(motifs))
	#group features by frame and locatiion
	groups = group_motifs(motifs, max_gap=1500)
	if verbose:
		print "Got {} groups, extracting envelopes...".format(len(groups))

	pprs = []
	dbg_env = []
	
	while groups:
		if verbose:
			print "Got {} groups, extracting envelopes...".format(len(groups))
		#extract the sequence envelope around each group
		envelopes = [get_envelope(group, target, margin=1000) for group in groups]
		dbg_env += envelopes

		if verbose:
			print "Got {} envelopes, locating PPRs...".format(len(envelopes))
		#locate the PPR within each envelope
		for envelope in envelopes:
			ppr = locate_ppr(envelope)
			if ppr:
				pprs.append(ppr)

		#look for overlapping pprs
		groups = remove_overlaps(pprs)
		ol = len(groups)
		if verbose:
			print "{} conflicts".format(ol)
		groups += remove_overgrown(pprs, 500)
		if verbose:
			print "{} overgrown PPRs".format(len(groups) - ol)
	
	pprs = [add_source(p, target) for p in pprs]
	
	if verbose:
		print "Got {} PPRs, cleaning...".format(len(pprs))
	#clean the gaps between features
	pprs = [clean_ends(fill_gaps(ppr)) for ppr in pprs]

	#annotate the tail region and classify each PPR
	classify.classify(pprs)

	#predict each PPR's target
	targetp.targetp(pprs, annotation='localization')

	#filter the desired location
	if localization:
		pprs = [p for p in pprs if p.annotations['localization'] == localization]

	#return a list of nicely presented PPRs
	return pprs

def debug_ends(ppr):
	for m in ppr.features:
		if m.location.start < 0:
			print "\tPPR motif begins before protein does!"

def group_motifs(motifs, max_gap):
	"""Given a list of PPR motifs, group them by strand then by location such
	that the maximum gap within a group is max_gap"""
	def distance(a,b):
		"""return the minimum distance between to seqfeatures"""
		return min(	abs(a.location.start - b.location.end),
								abs(b.location.start - a.location.end))

	strands = {1: [], -1: [],}
	for motif in motifs:
		if motif.qualifiers['frame'] > 0:
			strands[1].append(motif)
		else:
			strands[-1].append(motif)

	groups = []
	#for each frame
	for strand in strands.itervalues():
		#sort by start location
		strand.sort(key=lambda m: m.location.start)
		#for each motif
		current_group = []
		for motif in strand:
			#if the group is empty, initiate it
			if not current_group:
				current_group.append(motif)
				continue

			#if we're within the distance limit, add to the group
			if distance(motif, current_group[-1]) < max_gap:
				current_group.append(motif)
			#otherwise start a new one
			else:
				groups.append(current_group)
				current_group = [motif,]

		if current_group:
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
	strand = group[0].location.strand

	seq = target.seq[start:end]
	if strand < 0:
		seq = seq.reverse_complement()
	seq.alphabet=generic_dna

	return SeqRecord(seq, id='env', name='PPR_env', description='PPR_envelope',
			annotations = {'src_from': start, 'src_to': end, 
				'src_strand': strand})

def locate_ppr(envelope):
	"""Find and annotate the protein within"""
	#find all the PPR motifs
	search = HMMER.hmmsearch(hmm = models[3], targets = envelope)
	motifs = search.getFeatures(envelope)

	#A ppr must contain 2 or more PPR motifs
	if len(motifs) < 2:
		return None
	
	#order the motifs
	motifs.sort(key=lambda m: m.location.start)
	known_start = True
	known_stop = True
	#find start codon
	start = motifs[0].location.start
	while start > 0 and str(envelope.seq[start:start+3]).lower() != "atg":
		start -= 3
	if start < 0:
		known_start = False
		start = 0

	#find stop codon
	stop = motifs[-1].location.end
	while stop < len(envelope) and (
		str(envelope.seq[stop:stop+3]).lower() not in ["tag", "tga", "taa"]):
		stop += 3
	if stop > len(envelope):
		known_stop = False
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

	annotations = {
				"src_from"	: src_from,
				"src_to"		: src_to,
				"src_strand": envelope.annotations['src_strand'],
	}
	if not known_stop:
		annotations['no_stop'] = True
	if not known_start:
		annotations['no_start'] = True

	#return a record
	return SeqRecord(envelope.seq[start:stop],
			features = motifs,
			annotations = annotations)

def remove_overlaps(pprs):
	"""Remove any overlapping PPRs from pprs and return their groups"""
	def range_check(a,b,c):
		"""SUE ME!!!!"""
		return c >= a and c <= b
	
	def overlaps(a,b):
		"""return true if a overlaps with b on the source"""
		if a.annotations['src_strand'] != b.annotations['src_strand']:
			return False

		a_from = a.annotations['src_from']
		a_to   = a.annotations['src_to']
		b_from = b.annotations['src_from']
		b_to   = b.annotations['src_to']
		return (range_check(a_from, a_to, b_from) or
						range_check(a_from, a_to, b_to) or
						range_check(b_from, b_to, a_from))

	pprs.sort(key=lambda p: p.annotations['src_from'])

	conflicts = []
	current_conflict = []
	last_ppr = None
	for ppr in pprs:
		if current_conflict:
			if sum([overlaps(ppr, c) for c in current_conflict]) > 0:
				current_conflict.append(ppr)
			else:
				conflicts.append(current_conflict)
				current_conflict = None
				last_ppr = ppr
		elif last_ppr and overlaps(ppr, last_ppr):
			current_conflict = [last_ppr, ppr]
			last_ppr = None
		else:
			last_ppr = ppr

	#remove all conflicting PPRs from pprs
	for conflict in conflicts:
		for ppr in conflict:
			pprs.remove(ppr)

	#return groups from each PPR
	groups = []
	for conflict in conflicts:
		groups.append([SeqFeature(FeatureLocation(
			min([ppr.annotations['src_from'] for ppr in conflict]),
			max([ppr.annotations['src_to'] for ppr in conflict]),
			conflict[0].annotations['src_strand'])),]
		)

	return groups

def remove_overgrown(pprs, margin):
	"""Remove PPRs which turned out to be too large for the envelope"""
	overgrown = []

	#find pprs which are missing start or stop codons due to not having enough
	#space
	for p in pprs:
		if (p.annotations.has_key('no_start') and 
				p.features[0].location.start < margin):
			overgrown.append(p)
		elif (p.annotations.has_key('no_stop') and
				(len(p) - p.features[-1].location.end) < margin):
			overgrown.append(p)
	
	#remove the from pprs
	for p in overgrown:
		pprs.remove(p)
	
	#return groups from each PPR
	groups = []
	for ppr in overgrown:
		groups.append([SeqFeature(FeatureLocation(
			ppr.annotations['src_from'],
			ppr.annotations['src_to'],
			ppr.annotations['src_strand'])),]
		)

	return groups

def add_source(ppr, source):
	"""Add source metadata to the ppr"""
	ppr.annotations.update({
		'src_id': source.id,
		'src_name': source.name,
		'src_desc': source.description,
		})
	ppr.id = "PPR"
	ppr.name = "PPR_Protein"
	ppr.description = "PPR Protein"
	return ppr

def fill_gaps(ppr):
	"""Look for reluctant motifs in suspiciously sized gaps"""

	gaps = find_gaps(ppr, mingap=3*30, maxgap=3*40)	
	#check if there could be a motif at the start
	if ppr.features[0].location.start > 3*30:
		gaps.append(FeatureLocation(0, ppr.features[0].location.start, 1))

	for g in gaps:
		record = ppr[g.start-15:g.end+15]
		record.seq = record.seq.translate()

		#increase the reporting thresholds
		search = HMMER.hmmsearch(hmm = models[3], targets = record, 
				F1=0.5,F2=0.5,F3=0.5)
		
		if search.matches:
			motif = search.getFeatures(record)[0]
			offset = g.start-15
			motif.location = FeatureLocation(
					offset + 3*motif.location.start,
					offset + 3*motif.location.end,
					strand = 1)
			ppr.features.append(motif)
	ppr.features.sort(key=lambda p: p.location.start)
	return ppr

def clean_ends(ppr):
	"""Clean the length of the PPR motifs such that they start at the correct
	location and so that small gaps are eliminated

	adds the qualifier type as P, L or S
	"""
	gaps = find_gaps(ppr, mingap=None, maxgap=10)
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

def find_gaps(ppr, mingap=30, maxgap=None, skip_introns=True):
	"""Find all the gaps between PPR motifs which are gte maxgap"""
	loc = []
	feats = sorted(ppr.features, key = lambda(p): int(p.location.start))
	for a,b in pairwise(feats):
		#ignore any gaps if a and b aren't in the same frame
		if a.qualifiers['frame'] != b.qualifiers['frame']:
			continue

		#if the size is within the range
		l = int(b.location.start) - int(a.location.end)
		if l >= (mingap or -float('inf')) and l <= (maxgap or float('inf')):
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
		

