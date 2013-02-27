"""Useful utilities"""

import os, os.path, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyHMMER import HMMER, hmmfile

HMMDir = os.path.join(os.path.dirname(__file__), 'HMMs/')
TargetDir = os.path.join(os.path.dirname(__file__), 'Genomes/')
TestDir = os.path.join(os.path.dirname(__file__), 'Test_Proteins/')

def getLabelledFeatures(seq, search='pentatricopeptide', feat_type=None):
	"""Return all the features in seq which contain desc and are of type type"""
	if isinstance(seq, SeqRecord):
		lseq = [seq,]
	else:
		lseq = seq

	ret = []
	
	for i,s in enumerate(lseq):
		r = []
		for feat in s.features:
			d = str(feat.qualifiers.values())
			if search and not search.lower() in d:
				continue
			if feat_type and not feat_type.lower() == feat.type.lower():
				continue
			r.append(feat)
		ret.append(r)
	
	if isinstance(seq, SeqRecord):
		return ret[0]
	return ret

def loadmodels():
	"""Load all the models from a file"""
	modelfiles = [ os.path.join(HMMDir, "PPR_{}.hmm".format(i)) for i in
			range(4)]
	return [hmmfile.read(f)[0] for f in modelfiles]

def gettargetnames():
	"""get the filenames of all the available target genomes"""
	return [os.path.abspath(os.path.join(TargetDir, path)) for 
			path in os.listdir(TargetDir)]

def loadnuclear(f):
	"""Load the sequences of all the nuclear genomes found in f"""
	gen = _open_seq(f)	

	ignore = ['chloroplast', 'mitochondria',]
	def shouldignore(rec):
		for i in ignore:
			if (rec.id.lower().find(i) >= 0 or rec.name.lower().find(i) >= 0 or
					rec.description.lower().find(i) >= 0):
				return True
		return False

	ret = []
	for rec in gen:
		if not shouldignore(rec):
			ret.append(rec)

	return ret

def loadchloroplast(f):
	"""Load the chloroplast sequence from f"""
	c = 'chloroplast'
	for rec in _open_seq(f):
		if (rec.name.find(c) >=0 or rec.id.find(c) >= 0):
			return rec

def _open_seq(f):
	if os.path.splitext(f)[1].lower() in ['.gb', '.gbk', '.genbank', '.gen',]:
		return SeqIO.parse(f, 'genbank')
	elif os.path.splitext(f)[1].lower() in ['.fas', '.fasta',]:
		return SeqIO.parse(f, 'fasta')
	else:
		raise ValueError('Could not detect file type for \'{}\''.format(f))

def get_tail_consensus():
	return (SeqIO.read(os.path.join(HMMDir, 'E.fasta'), 'fasta'),
		SeqIO.read(os.path.join(HMMDir, 'E+.fasta'), 'fasta'),
		SeqIO.read(os.path.join(HMMDir, 'DYW.fasta'), 'fasta'))

def get_tail_models():
	return (hmmfile.read(os.path.join(HMMDir, 'E.hmm'))[0],
		hmmfile.read(os.path.join(HMMDir, 'E+.hmm'))[0],
		hmmfile.read(os.path.join(HMMDir, 'DYW.hmm'))[0])

def load_test():
	s = SeqIO.read(os.path.join(TestDir, 'PPR10.gb'), 'genbank')
	s.features = []
	return s

import itertools
def pairwise(iterable):
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = itertools.tee(iterable)
	next(b, None)
	return itertools.izip(a, b)
