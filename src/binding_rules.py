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


def build_yagi(ppr):
	NSRs = get_yagi_NSRs(ppr)
	PSSM = []
	for nsr in NSRs:
		PSSM.append(tuple((log(p) for p in NSR_to_row(nsr))))
	return PSSM

CODES['yagi'] = build_yagi

def get_yagi_NSRs(ppr):
	"""Return a list of 3-tuples referring to the NSRs"""
	o = 3 * 1
	NSR = []
	for a,b in pairwise(ppr.features):
		s1 = ppr.seq[a.location.start+o:o+a.location.start + 1 + 3 * 4].translate()
		s2 = ppr.seq[b.location.start - 3*2+o:o+b.location.start].translate()
		NSR.append((s1[0], s1[3], s2[0],))
	return NSR

def NSR_to_row(NSR):
	s = "".join(NSR).upper()
	if YAGI_3.has_key(s):
		return YAGI_3[s]
	s = '*' + s[1:]
	if YAGI_2.has_key(s):
		return YAGI_2[s]
	s = '*' + s[1] + '*'
	if YAGI_1.has_key(s):
		return YAGI_1[s]
	return [0.25, 0.25, 0.25, 0.25]

YAGI_3 = {
	'VND': [0.08, 0.31, 0.03, 0.58],
	'FND': [0.18, 0.20, 0.10, 0.52],
	'EGD': [0.07, 0.05, 0.69, 0.19],
	'VNN': [0.21, 0.61, 0.05, 0.13],
	'FSN': [0.63, 0.14, 0.09, 0.14],
	'VTN': [0.49, 0.18, 0.23, 0.10],
	'ITD': [0.08, 0.20, 0.47, 0.25],
	'INN': [0.08, 0.50, 0.05, 0.36],
	'YND': [0.09, 0.15, 0.14, 0.62],
	'TSD': [0.14, 0.18, 0.14, 0.54],
	'VSN': [0.53, 0.12, 0.06, 0.29],
	'FTN': [0.75, 0.07, 0.06, 0.12],
	'VNS': [0.09, 0.62, 0.06, 0.23],
	'FPD': [0.10, 0.08, 0.06, 0.77],
	'YPD': [0.10, 0.08, 0.06, 0.77],
	'VSD': [0.41, 0.08, 0.27, 0.24],
	'VTD': [0.10, 0.08, 0.70, 0.13],
	'VVE': [0.31, 0.18, 0.06, 0.45],
	'EGN': [0.52, 0.08, 0.06, 0.34],
	'FNN': [0.45, 0.08, 0.34, 0.13],
	'YTN': [0.10, 0.08, 0.27, 0.55],
	'FNW': [0.17, 0.29, 0.27, 0.27],
	'INW': [0.10, 0.29, 0.06, 0.55],
	'VNW': [0.20, 0.39, 0.06, 0.34],}

YAGI_2 = {
	'*ND': [0.10, 0.18, 0.11, 0.60],
	'*NN': [0.17, 0.47, 0.16, 0.21],
	'*TN': [0.57, 0.11, 0.15, 0.17],
	'*SD': [0.22, 0.16, 0.32, 0.29],
	'*TD': [0.16, 0.10, 0.59, 0.15],
	'*SN': [0.58, 0.10, 0.05, 0.26],
	'*GD': [0.10, 0.05, 0.60, 0.26],
	'*NW': [0.15, 0.34, 0.14, 0.36],
	'*PD': [0.07, 0.09, 0.10, 0.74],
	'*VN': [0.29, 0.22, 0.15, 0.34],
	'*NS': [0.10, 0.52, 0.04, 0.34],
	'*NT': [0.16, 0.50, 0.18, 0.16],
	'*GN': [0.55, 0.06, 0.05, 0.34],
	'*VG': [0.22, 0.48, 0.05, 0.25],
	'*AD': [0.14, 0.07, 0.56, 0.23],
	'*VD': [0.25, 0.24, 0.22, 0.29],
	'*VE': [0.25, 0.15, 0.06, 0.54],
	'*CD': [0.31, 0.15, 0.34, 0.20],
	'*AN': [0.24, 0.22, 0.06, 0.48],
	'*IN': [0.20, 0.29, 0.06, 0.45],
	'*LS': [0.17, 0.36, 0.06, 0.41],
	'*SS': [0.10, 0.29, 0.13, 0.48],
	'*VS': [0.31, 0.29, 0.06, 0.34],
	'*GT': [0.31, 0.29, 0.06, 0.34],
	'*TW': [0.41, 0.08, 0.17, 0.34],}

YAGI_1 = {
	'*N*': [0.12, 0.34, 0.12, 0.42],
	'*T*': [0.47, 0.09, 0.28, 0.16],
	'*S*': [0.36, 0.13, 0.21, 0.30],
	'*V*': [0.28, 0.22, 0.10, 0.39],
	'*G*': [0.31, 0.11, 0.30, 0.29],
	'*P*': [0.10, 0.06, 0.11, 0.72],
	'*A*': [0.27, 0.20, 0.18, 0.35],
	'*L*': [0.18, 0.27, 0.04, 0.52],
	'*I*': [0.20, 0.39, 0.10, 0.30],
	'*C*': [0.55, 0.10, 0.21, 0.14],
	'*M*': [0.10, 0.15, 0.06, 0.69],}


from Bio import SeqIO, Alphabet
import extract

def test_yagi_NSRs():
	"""Test yagi against PpPPR_71"""

	for ppr in SeqIO.parse("Test_Proteins/Known_PPRs.gb", 'genbank'):
		print "{}: {}".format(ppr.name,	NSR_to_string(get_yagi_NSRs(ppr)))

def print_motif(ppr, number):
	f = ppr.features[number-1]
	print "Motif {}: {}|{}|".format(number, 
			ppr.seq[f.location.start-6:f.location.start+1].translate(),
			f.extract(ppr).seq.translate())

def NSR_to_string(NSR):
	s = ""
	for n in NSR:
		s += ''.join(n).upper()
		if n != NSR[-1]:
			s += '-'
	return s


