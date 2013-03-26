"""Useful utilities"""

import os, os.path, re
from Bio import SeqIO, Alphabet
from Bio.SeqRecord import SeqRecord
from pyHMMER import HMMER, hmmfile

HMMDir = os.path.join(os.path.dirname(__file__), 'HMMs/')
TargetDir = os.path.join(os.path.dirname(__file__), 'Genomes/')
TestDir = os.path.join(os.path.dirname(__file__), 'Test_Proteins/')
TestData = os.path.join(os.path.dirname(__file__), 'Test_Data/')

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

	for rec in gen:
		if not shouldignore(rec):
			yield rec

def loadchloroplast(f):
	"""Load the chloroplast sequence from f"""
	c = 'chloroplast'
	for rec in _open_seq(f):
		if (rec.name.find(c) >=0 or 
				rec.id.find(c) >= 0 or 
				rec.description.find(c) >= 0):
			return rec

def _open_seq(f):
	if os.path.splitext(f)[1].lower() in ['.gb', '.gbk', '.genbank', '.gen',]:
		return SeqIO.parse(f, 'genbank', alphabet=Alphabet.generic_dna)
	elif os.path.splitext(f)[1].lower() in ['.fas', '.fasta',]:
		return SeqIO.parse(f, 'fasta', alphabet=Alphabet.generic_dna)
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

def load_test_targets():
	return list(SeqIO.parse(os.path.join(TestData, 'PPR10_targets.fas'), 'fasta')) 

import itertools
def pairwise(iterable):
	"s -> (s0,s1), (s1,s2), (s2, s3), ..."
	a, b = itertools.tee(iterable)
	next(b, None)
	return itertools.izip(a, b)

def write_data(header, data, filename):
	"""
		Write data to file in PGF plots friendly formatting
		header = column headers
		data = list of tuples (rows)
		filename = file to write to
	"""
	#Data and header must be present
	if not data:
		raise ValueError("No data given")
	if not header:
		raise ValueError("Header empty")

	#each data item and the header must have the same number of items
	fields = len(header)
	for i,d in enumerate(data, 1):
		if len(d) != fields:
			raise ValueError("Data row {} has {} too {} fields".format(
				i, abs(len(d)-fields), "many" if len(d)>fields else "few"))
	
	#each column must be entirely the same type
	types = tuple((type(d) for d in data[0]))
	for i,d in enumerate(data,1):
		for j,item in enumerate(d):
			if not isinstance(item, types[j]):
				raise TypeError("\'{}\' in row {} is of type {} not {}".format(
					header[j], i, type(item), types[j]))

	formats = {
			str: "{{:<{}s}}",
			int: "{{:>{}d}}",
			float: "{{:>{}.4f}}",
	}
	hformats = {
			str: "{{:<{}s}}",
			int: "{{:>{}s}}",
			float: "{{:>{}s}}",
	}

	#Calculate the maximum width of each
	widths = []
	for i in range(len(header)):
		w = len(formats[str].format(header[i]))
		widths.append(max(w, 
			*[len(formats[types[i]].format(1).format(d[i])) for d in data]))

	try:
		header_fmt = (" "
				.join([hformats[t].format(w) for t,w in zip(types,widths)]) + '\n')
		row_fmt = (" "
				.join([formats[t].format(w) for t,w in zip(types,widths)]) + '\n')
	except KeyError, e:
		raise ValueError("Unsupported type {}".format(e.message))

	f = open(filename, 'w')
	f.write(header_fmt.format(*header))
	for d in data:
		f.write(row_fmt.format(*d))

	f.close()

	return len(data)
	


