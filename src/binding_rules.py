"""Given a PPR SeqRecord, build a PSSM and return it"""

from utils import pairwise
from math import log, exp, sqrt, pi

equal = (1, 1, 1, 1)

CODES={}

def build_PSSM(ppr, coding='barkan', **kwargs):
	"""Build a log-odds PSSM for the ppr's target Vs background"""

	if CODES.has_key(coding.lower()):
		return CODES.get(coding)(ppr, **kwargs)
	else:
		raise ValueError("Unknown coding \'{}\', only {} are accepted".format(
				coding, ', '.join(CODES.iterkeys())))
	
def build_barkan(ppr, background=None):
	"""Build a PSSM based on Barkan coding"""

	if background == None:
		background = (1,1,1,1,)
	code = get_code(ppr)

	#print "Building model for:\n\t{}".format(
	#		string_code(code).replace('\n','\n\t'))

	PSSM = []

	for i,(a,b) in enumerate(pairwise(code)):
		s = a[1] + b[0]
		if a[2] == 'P' and s in Ptype:
			emit = Ptype[s]
		elif a[2] == 'S' and s in Stype:
			emit = Stype[s]
		else:
			emit = equal

		if sum(emit) == 0:
			emit = equal


		#convert to log odds
		tot = sum(emit)
		emit = tuple(log(float(i) / float(tot*b)) for i,b in zip(emit,background))
		#print "{:2}: \"{}\" [{}] -> {}".format(i,s,a[2],emit)
		PSSM.append(emit)

	return PSSM

CODES['barkan'] = build_barkan

def get_code(ppr):
	"""Return the amino acids at the 1 and 6 positions

	-  ppr: a SeqRecord with annotated PPR repeat domains 

	-  return: [(1_1, 1_6, type,), (2_1, 2_6, type), ...]
	"""
	ret = []
	for motif in ppr.features:
		start = int(motif.location.start)
		seq = ppr.seq[start:start+18].translate()
		type = motif.qualifiers['type']
		if not isinstance(type, basestring):
			type = type[0]
		ret.append( (seq[0],seq[5], type,) )
	return ret

def string_code(code):
	"""get a human readable representation"""
	return "6 : {}\n1 :{}\nT : {}".format(
			"".join([x[1] for x in code]),
			"".join([x[0] for x in code]),
			"".join([x[2] for x in code]))


Ptype = {
		#		  A   C   G   T
	'TD': ( 2,  0, 23,  1, ),
	'TN': (23,  0,  2,  0,),
	'ND': ( 7, 22,  7, 68,),
	'NN': ( 6, 46,  3, 31,),
	'SN': (10,  1,  0,  1,),
	'NS': ( 1, 12,  0,  3,),
	'TT': ( 0,  0,  3,  0,),
	'SS': ( 5,  1,  0,  0,),
	'SR': ( 0,  0,  2,  0,),
	'SC': ( 2,  0,  0,  0,),
	'RD': ( 0,  0,  0,  3,),
	'TS': ( 3,  0,  0,  1,),
	'NR': ( 0,  2,  0,  0,),
	'AD': ( 0,  0,  2,  1,),
	'GD': ( 0,  1,  2,  0,),
	'HM': ( 0,  0,  0,  2,),
	'HS': ( 0,  0,  0,  2,),
	'NC': ( 0,  5,  0,  3,),
	'FG': ( 1,  0,  0,  0,),
	'GH': ( 1,  0,  0,  0,),
	'GS': ( 1,  0,  0,  0,),
	'SG': ( 1,  0,  0,  0,),
	'SL': ( 1,  0,  0,  0,),
	'TP': ( 1,  0,  0,  0,),
	'HV': ( 2,  0,  0,  1,),
	'GN': ( 0,  1,  0,  0,),
	'NT': ( 1,  8,  0,  8,),
	'NG': ( 1,  2,  0,  5,),
	'AN': ( 2,  2,  0,  0,),
	'MD': ( 0,  1,  1,  0,),
	'NE': ( 0,  0,  1,  1,),
	'TR': ( 1,  0,  1,  0,),
	'SD': ( 4,  3,  3,  1,),
	'CS': ( 0,  1,  0,  2,),
	'ID': ( 1,  0,  1,  1,),
}

Stype = {
		#		  A   C   G   T
	'TD': ( 9,  2, 27,  4,),
	'SD': ( 1,  1, 13,  1,),
	'SN': (20,  2,  1,  5,),
	'TN': (18,  1,  2,  3,),
	'ND': ( 5, 13,  5, 44,),
	'CN': ( 3,  0,  0,  0,),
	'NT': ( 0, 12,  3,  4,),
	'NN': ( 5,  8,  9,  5,),
	'GS': ( 2,  0,  0,  0,),
	'AR': ( 0,  2,  0,  0,),
	'TH': ( 1,  0,  2,  0,),
	'CD': ( 0,  0,  1,  0,),
	'GD': ( 0,  0,  1,  0,),
	'PC': ( 0,  0,  1,  0,),
	'AH': ( 1,  0,  0,  0,),
	'SE': ( 1,  0,  0,  0,),
	'TE': ( 1,  0,  0,  0,),
	'TK': ( 1,  0,  0,  0,),
	'SK': ( 2,  0,  1,  0,),
	'TL': ( 2,  0,  0,  1,),
	'LD': ( 0,  1,  0,  0,),
	'NE': ( 0,  1,  0,  0,),
	'AT': ( 0,  0,  0,  1,),
	'FN': ( 0,  0,  0,  1,),
	'IN': ( 0,  0,  0,  1,),
	'PN': ( 0,  0,  0,  1,),
	'SS': ( 0,  0,  0,  1,),
	'TQ': ( 0,  0,  0,  1,),
	'VS': ( 0,  0,  0,  1,),
	'KD': ( 0,  2,  1,  0,),
	'NH': ( 0,  2,  0,  1,),
	'AD': ( 0,  0,  1,  1,),
	'TS': ( 1,  3,  1,  1,),
	'AN': ( 1,  0,  0,  1,),
	'TT': ( 1,  1,  1,  0,),
	'NS': ( 2,  4,  1,  5,),
}

#add in a prior
for k,v in Ptype.iteritems():
	Ptype[k] = tuple(i + 1 for i in v)
for k,v in Stype.iteritems():
	Stype[k] = tuple(i + 1 for i in v)
