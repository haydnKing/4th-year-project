"""This file tests the code and pyHMMER against the annotated version of the
	arabidopsis genome found in EBI"""

from Bio import SeqIO
import extract
import os.path

DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(DIR, 'Test Data/TAIR10_PPRs.gb')
ARA_FILE = os.path.join(DIR, 'Genomes/Arabidopsis_TAIR10.fas')

def get_known():
	l = list(SeqIO.parse(DATA_FILE, 'genbank'))
	return [f.features for f in l]

def print_numbers_found(found, known):
	fmt = "{: <20} {: <10} {: <10}"
	print fmt.format('Chromosome', 'Found', 'Known')

	for i,(f,k) in enumerate(zip(found,known)):
		print fmt.format(i+1, len(f), len(k))
	print fmt.format('Sigma', sum([len(f) for f in found]), sum([len(k) for k in
		known]))

def find_match(f, known):
	"""Search in known for the PPR found"""
	p = f.annotations['sourceloc']
	p = int((int(p.start) + int(p.end))/2)
	for k in known:
		if p in k:
			return k
	return None

def compare_chromosomes(found_pprs, known):
	found = []
	missed = []
	duplicate = []
	false = []

	for f in found_pprs:
		k = find_match(f, known)
		if not k:
			false.append(f)
		else:
			if not hasattr(k, 'found'):
				k.hasattr = True
				found.append(f)
			else:
				duplicate.append(f)
	
	for k in known:
		if not hasattr(k, 'found'):
			missed.append(k)
		else: 
			delattr(k, 'found')

	print ("\tfound: {}\n\tmissed: {}\n\tduplicate: {}\n\tfalse: {}"
		.format(len(found), len(missed), len(duplicate), len(false)))

def get_found():
	ara = list(SeqIO.parse(ARA_FILE, 'fasta'))
	pprs = extract.simple_extract(ara)

	found = [[],[],[],[],[]]
	for ppr in pprs:
		try:
			found[ara.index(ppr.annotations['source'])].append(ppr)
		except IndexError:
			pass

	return found

def compare():

	found = get_found()
	known = get_known()

	print_numbers_found(found, known)
	for i,(f,k) in enumerate(zip(found,known)):
		print "Chromosome {}".format(i+1)
		compare_chromosomes(f,k)

if __name__ == "__main__":
	compare()


