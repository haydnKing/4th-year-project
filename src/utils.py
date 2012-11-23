"""Useful utilities"""

import os, os.path, re
from Bio import SeqIO

HMMDir = os.path.join(os.path.dirname(__file__), 'HMMs/')
TargetDir = os.path.join(os.path.dirname(__file__), 'Genomes/')
TestDir = os.path.join(os.path.dirname(__file__), 'Test Proteins/')

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
		if (rec.name.find(c) >=0 or rec.id.find(c)):
			return rec

def _open_seq(f):
	if os.path.splitext(f)[1].lower() in ['.gb', '.gbk', '.genbank', '.gen',]:
		return SeqIO.parse(f, 'genbank')
	elif os.path.splitext(f)[1].lower() in ['.fas', '.fasta',]:
		return SeqIO.parse(f, 'fasta')
	else:
		raise ValueError('Could not detect file type for \'{}\''.format(f))
