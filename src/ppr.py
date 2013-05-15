"""PPR utils"""

import PSSM, os, json, re
from Bio import SeqIO

class PPR:
	"""Represent a PPR"""
	
	def __init__(self, seq_record):
		self.seq_record = seq_record
		self.pssm = PSSM.PSSM_build(self.seq_record)

	def distance(self, rhs, method='MSE'):
		"""Calculate the distance between myself and another PPR using method:
			method =	'MSE': Mean Squared Error
								'KL': KL divergance
		"""
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

	def distance(self, ppr, method=None):
		if method:
			return [ppr.distance(p, method)[0] for p in self.pprs]
		return [ppr.distance(p)[0] for p in self.pprs]

	def closest(self, ppr, method=None):
		r = self.distance(ppr, method)
		return self.pprs[r.index(min(r))]

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


