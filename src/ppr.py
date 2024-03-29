"""PPR utils"""

import PSSM, os, json, re, extract
from Bio import SeqIO, Entrez
Entrez.email="hjk38@cam.ac.uk"

class PPR:
	"""Represent a PPR"""
	
	def __init__(self, seq_record):
		self.seq_record = seq_record
		self.pssm = PSSM.build(self.seq_record)

	def distance(self, rhs, method='MSE'):
		"""Calculate the distance between myself and another PPR using method:
			method =	'MSE': Mean Squared Error
								'KL': KL divergance
		"""
		if method == None:
			method = "MSE"
		return PSSM.distance(self.pssm, rhs.pssm, method)

	def localization(self):
		"""Return the TargetP localization prediction"""
		return self.seq_record.annotations.get('localization', "*").upper()

	def type(self):
		"""Return the PPR's type (P,L or S)"""
		return self.seq_record.annotations.get('type', "").upper()

	def family(self):
		return self.seq_record.annotations.get('ppr_family', "").upper()

	def __len__(self):
		return len(self.seq_record.features)

	def __repr__(self):
		return "<PPR type={}, localization={}, length={}>".format(
				self.type(), self.localization(), len(self))

	def __str__(self):
		return self.__repr__()


class Genome:
	"""Represent a Genome"""
	def __init__(self, name):
		self.name = name
		self.pprs = load_pprs(name)

	def __str__(self):
		return "Genome {} ({} pprs)".format(self.name, len(self.pprs))

	def __len__(self):
		return len(self.pprs)

	def distance(self, ppr, method=None):
		if method:
			return [ppr.distance(p, method)[0] for p in self.pprs]
		return [ppr.distance(p)[0] for p in self.pprs]

	def closest(self, ppr, method=None):
		r = self.distance(ppr, method)
		if r:
			return self.pprs[r.index(min(r))]
		else:
			return None

def get_genomes():
	genomes = [f for f in os.listdir('output/PPRs/') if re.match(".+\.gb$", f)]
	return sorted([os.path.splitext(f)[0] for f in genomes])

def load_pprs(genome):
	return [PPR(p) for p in load_records(genome)]

def load_records(genome):
	path = "output/PPRs/{}.gb".format(genome)
	annot = "output/PPRs/{}.json".format(genome)
	pprs = list(SeqIO.parse(path, "genbank"))
	f = open(annot, "r")
	for a,p in zip(json.loads(f.read()), pprs):
		p.annotations = a
	f.close()
	return pprs		

def load_genomes():
	return [Genome(name) for name in get_genomes()]

def write_pprs(pprs, genome):
	out = "output/PPRs/{}.gb".format(genome)
	annot = "output/PPRs/{}.json".format(genome)
	SeqIO.write(pprs, out, "genbank")
	annots = []
	for p in pprs:
		annots.append(p.annotations)
	f = open(annot, "w")
	f.write(json.dumps(annots))
	f.close()

def save_known_pprs():
	names = ["PDE247","CRR21","CRR22","CRR28","CRR4","LPA66","otp80","otp81",
			"otp82","otp84","otp85","otp86",]
	pprs = []
	for n in names:
		print "Search \'{}\'".format(n)
		hnd = Entrez.esearch(db='nuccore', 
				term='({}[Title]) AND Arabidopsis thaliana[Organism]'.format(n))
		ret = Entrez.read(hnd)
		if len(ret['IdList']) != 1:
			print "failed to find \'{}\'".format(n)
		else:
			p = extract.from_entrez(ret['IdList'][0])[0]
			p.name = n
			pprs.append(p)
	
	print "Found {} PPRs".format(len(pprs))
	SeqIO.write(pprs, "Test_Proteins/Known_PPRs.gb", "genbank")

