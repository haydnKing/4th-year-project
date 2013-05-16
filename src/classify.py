"""Classify PPR proteins"""

import extract, utils, os.path
from pyHMMER import HMMER, hmmfile
from Bio import SeqIO, Alphabet
from Bio.SeqRecord import SeqRecord

from StringIO import StringIO

def update_models():
	"""Recalculate the HMM models"""
	print "Extracting C-Termini..."
	ct = extract.get_c_terminus(extract.extract(localization=None))
	print "Done. Got {} tails".format(len(ct))
	
	(E, Ep, DYW) = utils.get_tail_consensus()

	print "E"
	j = HMMER.jackhmmer(E,ct)
	j.hmms[-1].name = 'E'
	hmmfile.write(j.hmms[-1], os.path.join(utils.HMMDir, 'E.hmm'))

	print "E+"
	j = HMMER.jackhmmer(Ep,ct)
	j.hmms[-1].name = 'E+'
	hmmfile.write(j.hmms[-1], os.path.join(utils.HMMDir, 'E+.hmm'))

	print "DYW"
	j = HMMER.jackhmmer(DYW,ct)
	j.hmms[-1].name = 'DYW'
	hmmfile.write(j.hmms[-1], os.path.join(utils.HMMDir, 'DYW.hmm'))

def classify(pprs, family_annot="ppr_family", tail_annot='ppr_tail'):
	"""Annotate each ppr with it's family type, (P,PLS,E,E+,DYW)"""
	ct = get_c_terminus(pprs)
	(E, Ep, DYW) = utils.get_tail_models()

	for i,c in enumerate(ct):
		if not isinstance(c.seq.alphabet, Alphabet.ProteinAlphabet):
			print "ct[{}]: {}".format(i, str(c.seq))

	h = HMMER.hmmsearch([E,Ep,DYW], ct)

	#annotate each tail
	h.annotate(ct)

	for ppr,tail in zip(pprs, ct):
		fmt = ''
		if tail.features:
			f = sorted(tail.features, key=lambda(ft): int(ft.location.start))
			fmt = ("-{.type}"*len(f)).format(*f)
			if fmt[-3:] == 'DYW':
				ppr.annotations[family_annot] = 'DYW'
			elif fmt.find('E+') >= 0:
				ppr.annotations[family_annot] = 'E+'
			elif fmt.find('E') >= 0:
				ppr.annotations[family_annot] = 'E'
			else:
				print "Unknown tail format \'{}\'".format(fmt)
				ppr.annotations[family_annot] = '??'

		else:
			l = len(ppr.features[0])
			for f in ppr.features:
				if len(f) != l:
					ppr.annotations[family_annot]='PLS'
					continue
			ppr.annotations[family_annot] = 'P'

		ppr.annotations[tail_annot] = fmt

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

	return ret

def print_hist(pprs):
	"""Display a histogram of how many pprs are in each class"""
	hist = {'P':0,'PLS':0,'E':0,'E+':0,'DYW':0}

	for p in pprs:
		try:
			c = p.annotations['ppr_family']
			hist[c] += 1
		except KeyError:
			pass
	
	for k,v in hist.iteritems():
		print "{:4s}: {}".format(k, '*'*v)

