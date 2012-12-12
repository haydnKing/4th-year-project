"""Produce an HMM predicting the binding specificity of an HMM"""

from pyHMMER.hmmfile import HMM



def get_1_6(ppr):
	"""Return the amino acids at the 1 and 6 positions

	-  ppr: a SeqRecord with annotated PPR repeat domains 

	-  return: [(1_1, 6_1,), (2_1, 6_1,), ...]
	"""
	ret = []
	for motif in ppr.features:
		start = int(motif.location.start)
		seq = ppr.seq[start:start+15].translate()
		ret.append( (seq[0],seq[5]) )
	return ret

def print_1_6(code):
	"""print out a human readable representation"""
	print "6 : {}".format("".join([x[1] for x in code]))
	print "1 :{}".format("".join([x[0] for x in code]))

def build_model(code):
	"""Build a model from the list of 1 and 6 amino acids"""
	pass

def get_binding(pprs):
	"""return a list of models where the nth model refers to the predicted
	binding target of the nth ppr in the list

	- pprs: a list of SeqRecords representing annotated PPR proteins

	- return: a list of pyHMMER.hmmfile.HMM objects
	"""
	pass

Ptype = {
	'TD': {'A': 2,  'C': 0 , 'G': 23, 'T': 1,},
	'TN': {'A': 23, 'C': 0 , 'G': 2 , 'T': 0 ,},
	'ND': {'A': 7 , 'C': 22, 'G': 7 , 'T': 68,},
	'NN': {'A': 6 , 'C': 46, 'G': 3 , 'T': 31,},
	'SN': {'A': 10, 'C': 1 , 'G': 0 , 'T': 1 ,},
	'NS': {'A': 1 , 'C': 12, 'G': 0 , 'T': 3 ,},
	'TT': {'A': 0 , 'C': 0 , 'G': 3 , 'T': 0 ,},
	'SS': {'A': 5 , 'C': 1 , 'G': 0 , 'T': 0 ,},
	'SR': {'A': 0 , 'C': 0 , 'G': 2 , 'T': 0 ,},
	'SC': {'A': 2 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'RD': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 3 ,},
	'TS': {'A': 3 , 'C': 0 , 'G': 0 , 'T': 1 ,},
	'NR': {'A': 0 , 'C': 2 , 'G': 0 , 'T': 0 ,},
	'AD': {'A': 0 , 'C': 0 , 'G': 2 , 'T': 1 ,},
	'GD': {'A': 0 , 'C': 1 , 'G': 2 , 'T': 0 ,},
	'HM': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 2 ,},
	'HS': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 2 ,},
	'NC': {'A': 0 , 'C': 5 , 'G': 0 , 'T': 3 ,},
	'FG': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'GH': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'GS': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'SG': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'SL': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'TP': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0 ,},
	'HV': {'A': 2 , 'C': 0 , 'G': 0 , 'T': 1 ,},
	'GN': {'A': 0 , 'C': 1 , 'G': 0 , 'T': 0 ,},
	'NT': {'A': 1 , 'C': 8 , 'G': 0 , 'T': 8 ,},
	'NG': {'A': 1 , 'C': 2 , 'G': 0 , 'T': 5 ,},
	'AN': {'A': 2 , 'C': 2 , 'G': 0 , 'T': 0 ,},
	'MD': {'A': 0 , 'C': 1 , 'G': 1 , 'T': 0 ,},
	'NE': {'A': 0 , 'C': 0 , 'G': 1 , 'T': 1 ,},
	'TR': {'A': 1 , 'C': 0 , 'G': 1 , 'T': 0 ,},
	'SD': {'A': 4 , 'C': 3 , 'G': 3 , 'T': 1 ,},
	'CS': {'A': 0 , 'C': 1 , 'G': 0 , 'T': 2 ,},
	'ID': {'A': 1 , 'C': 0 , 'G': 1 , 'T': 1 ,},
}

Stype = {
	'TD': {'A': 9 , 'C': 2 , 'G': 27, 'T': 4,},
	'SD': {'A': 1 , 'C': 1 , 'G': 13, 'T': 1,},
	'SN': {'A': 20, 'C': 2 , 'G': 1 , 'T': 5,},
	'TN': {'A': 18, 'C': 1 , 'G': 2 , 'T': 3,},
	'ND': {'A': 5 , 'C': 13, 'G': 5 , 'T': 44,},
	'CN': {'A': 3 , 'C': 0 , 'G': 0 , 'T': 0,},
	'NT': {'A': 0 , 'C': 12, 'G': 3 , 'T': 4,},
	'NN': {'A': 5 , 'C': 8 , 'G': 9 , 'T': 5,},
	'GS': {'A': 2 , 'C': 0 , 'G': 0 , 'T': 0,},
	'AR': {'A': 0 , 'C': 2 , 'G': 0 , 'T': 0,},
	'TH': {'A': 1 , 'C': 0 , 'G': 2 , 'T': 0,},
	'CD': {'A': 0 , 'C': 0 , 'G': 1 , 'T': 0,},
	'GD': {'A': 0 , 'C': 0 , 'G': 1 , 'T': 0,},
	'PC': {'A': 0 , 'C': 0 , 'G': 1 , 'T': 0,},
	'AH': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0,},
	'SE': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0,},
	'TE': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0,},
	'TK': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 0,},
	'SK': {'A': 2 , 'C': 0 , 'G': 1 , 'T': 0,},
	'TL': {'A': 2 , 'C': 0 , 'G': 0 , 'T': 1,},
	'LD': {'A': 0 , 'C': 1 , 'G': 0 , 'T': 0,},
	'NE': {'A': 0 , 'C': 1 , 'G': 0 , 'T': 0,},
	'AT': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'FN': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'IN': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'PN': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'SS': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'TQ': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'VS': {'A': 0 , 'C': 0 , 'G': 0 , 'T': 1,},
	'KD': {'A': 0 , 'C': 2 , 'G': 1 , 'T': 0,},
	'NH': {'A': 0 , 'C': 2 , 'G': 0 , 'T': 1,},
	'AD': {'A': 0 , 'C': 0 , 'G': 1 , 'T': 1,},
	'TS': {'A': 1 , 'C': 3 , 'G': 1 , 'T': 1,},
	'AN': {'A': 1 , 'C': 0 , 'G': 0 , 'T': 1,},
	'TT': {'A': 1 , 'C': 1 , 'G': 1 , 'T': 0,},
	'NS': {'A': 2 , 'C': 4 , 'G': 1 , 'T': 5,},
}

