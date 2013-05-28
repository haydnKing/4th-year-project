"""Test the binding rules on characterised PPRs"""

import binding, PSSM, utils, os.path, gc
from Bio import SeqIO, Alphabet
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import simplejson as json
from pyHMMER import HMMER

exact_codes = {
		'a': (10.0,  1.0,  1.0,  1.0,),
		'c': ( 1.0, 10.0,  1.0,  1.0,),
		'g': ( 1.0,  1.0, 10.0,  1.0,),
		't': ( 1.0,  1.0,  1.0, 10.0,),
}

def build_exact(footprint):
	"""Build a PSSM based on p.footprint"""
	PSSM = []
	for char in footprint.lower():
		PSSM.append(exact_codes[char])
	
	return PSSM

def find_exact(ppr, plastid):
	feats = []
	for e in ppr.exact:
		pfeats = binding.get_domains(e, plastid, 
				type="{}_exact".format(ppr.name), percentile=100.0, gaps=0)
		for feat in pfeats:
			if feat.qualifiers['odds'] >= 10.0 * len(e):
				feats.append(feat)
	return feats

def load_pprs():
	pprs = list(SeqIO.parse(os.path.join(utils.TestDir, "Known_PPRs.gb"), "gb",
							alphabet=Alphabet.generic_dna))
	
	#apply known binding sites
	f = open(os.path.join(utils.TestDir, "known_sites.json"))
	sites = json.load(f)
	for name,location in sites.iteritems():
		for p in pprs:
			if p.name.lower() == name.lower():
				p.footprints = [l.replace('u','t') for l in location]

	#build PSSMs
	for p in pprs:
		p.barkan = PSSM.build(p, 'barkan')
		p.yagi   = PSSM.build(p, 'yagi')
		p.exact  = [build_exact(footprint) for footprint in p.footprints]

	return pprs

def load_default_plastid():
	return SeqIO.read("Plastids/Arabidopsis_thaliana.gb", 'gb')

def annotate(pprs, plastid):
	for p in pprs:
		print "Annotate {}".format(p.name)
		print "  - yagi"
		feats =  binding.get_domains(p.yagi, plastid, 
				type="{}_yagi".format(p.name), gaps=2)
		
		print "  - barkan"
		feats += binding.get_domains(p.barkan, plastid, 
				type="{}_barkan".format(p.name), gaps=2)
		
		print "  - exact"
		feats += find_exact(p, plastid)

		plastid.features += feats

def save_known():
	"""Predict and save binding locations of all known PPR-RNAs in the
	arabidopsis chloroplast"""
	pprs = load_pprs()
	plastid = load_default_plastid()
	annotate(pprs, plastid)
	SeqIO.write(plastid, "output/ARA_annotated.gb", 'gb')

def load_plastids(exclude=[]):
	exclude_files = ["{}.gb".format(f.replace(' ','_')) for f in exclude]
	files = [f for f in os.listdir('Plastids') if (f not in exclude_files and
											not f.startswith('.'))]
	ret = []
	for f in files:
		seq = SeqIO.read("Plastids/{}".format(f), "gb",
				alphabet=Alphabet.generic_dna)
		seq.name = f[0:f.rfind('.')].replace('_',' ')
		ret.append(seq)

	return ret

def get_closest_gene(feature, genes):
	#filter by strand
	limit = [g for g in genes if g.location.strand == feature.location.strand]

	def dist(pos,f):
		ds = f.location.start - pos
		de = pos - f.location.end
		return min(ds,de)

	dists = [dist(feature.location.start, g) for g in limit]
	return limit[dists.index(min(dists))]

def sequence_similarity(s1, s2):
	r = 0
	for a,b in zip(str(s1).lower(),str(s2).lower()):
		if a == b:
			r += 1
	return float(r) / float(len(s1))

def find_homologs():
	"""Predict homologs of PPRs in other genomes based on footprints"""
	pprs = load_pprs()
	plastids = load_plastids(exclude=["Arabidopsis thaliana",])
	known_binding = SeqIO.read("output/ARA_annotated.gb", "gb")
	exact_features = [f for f in known_binding.features if 
																								"exact" in f.type.lower()]
	ara_genes = [f for f in known_binding.features if f.type.lower() == "gene"]
	ara_genes.sort(key=lambda g: g.location.start)

	print "Loaded {} pprs and {} plastids".format(len(pprs), len(plastids))

	for k,ppr in enumerate(pprs):
		print "Searching for homologs of \'{}\' ({}/{})".format(
				ppr.name,k+1,len(pprs))
		footprints = [f for f in exact_features if 
											f.type.lower() == "{}_exact".format(ppr.name.lower())]
		ppr.genes = [get_closest_gene(f, ara_genes) for f in footprints]

		print "\tFound {} original genes, {}".format(len(ppr.genes), 
				[g.qualifiers['gene'] for g in ppr.genes])

		ppr.potentialHomologs = {}

		for i,plastid in enumerate(plastids):

			if plastid.name != "Alsophila spinulosa":
				continue

			print "\t\tSearch {}/{}".format(i+1, len(plastids))

			#search for homologs of each gene
			homologs = []
			for gene in ppr.genes:
				g = SeqRecord(gene.extract(known_binding.seq).translate())
				search = HMMER.jackhmmer(g, plastid)
				print "{} -> {} homologs".format(gene.qualifiers['gene'],
						len(search.matches))
				homologs += search.getFeatures(type="{}_hl".format(gene.qualifiers['gene']))
			
			#extract the sequence surrounding each homolog
			for h in homologs:
				h.location = FeatureLocation(
						max(0,h.location.start - 500),
						min(len(plastid), h.location.end+500))
			homologs = [SeqRecord(h.extract(plastid.seq)) for h in homologs]

			#find exact or close to exact binding domains for each and add to the
			#list of potential homologs for the PPR
			ph = []
			for h in homologs:
				domains = []
				for exact in ppr.exact:
					try:
						domains += binding.get_domains(exact, h, percentile=100.0, gaps=0)
					except KeyError:
						continue
				if domains:
					domains.sort(key=lambda d: -d.qualifiers['odds'])
					seq = str(domains[0].extract(h).seq)
					similarity = max([sequence_similarity(original, seq) for 
																									original in ppr.footprints])
					print "  {} -> \'{}\'".format(h.type, seq)
					ph.append((similarity, seq))

			ph.sort(key=lambda p: -p[0])
			ppr.potentialHomologs[plastid.name] = ph
			
			#try and avoid running out of RAM
			gc.collect()

	for ppr in pprs:
		print "\'{}\' footprints = {}".format(ppr.name, ppr.footprints)
		print "potential homologs"
		for key,value in ppr.potentialHomologs.iteritems():
			print "{}: {}".format(key, value)

	return

	
	stats = []
	for plastid in plastids:
		length = 0
		similarity = 0.0
		for ppr in pprs:
			length += len(ppr.potentialHomologs[plastid.name])
			similarity += sum([p[0] for p in ppr.potentialHomologs[plastid.name]])
		
		try:
			stats.append({'name': plastid.name,
									'avg_similarity': similarity / float(length),
									'avg_homologs': length / len(pprs),})
		except ZeroDivisionError:
			stats.append({'name': plastid.name,
										'avg_similarity': 0.0,
										'avg_homologs': 0,})

	stats.sort(key=lambda s: -s['avg_similarity'])

	f = open("tmp", "w")

	for s in stats[0:50]:
		f.write("{name}, {avg_similarity}, {avg_homologs}\n".format(**s))
	f.close()





