"""Write out statistics from already processed genomes"""

import utils, os, os.path, ppr, classify

def length_hist(pprs):
	if not pprs:
		return [0,]
	lengths = [len(ppr) for ppr in pprs]
	hist = [0,] * (1+max(lengths))
	for l in lengths:
		hist[l] += 1
	return hist

def fmt_hist(hists):
	max_length = max([len(h) for h in hists]) - 1
	data = []
	for l in range(2, max_length+1):
		row = [l,]
		for h in hists:
			try:
				row.append(h[l])
			except IndexError:
				row.append(0)
		data.append(tuple(row))
	return data

def get_numbers(pprs):
	return (len(pprs),)

def get_numbers_hdr():
	return ("number",)

def get_localization(genome):
	data = []
	for g in genome:
		if len(g.pprs) < 50:
			continue
		total = float(len(g.pprs))
		loc_c = 100.0 * len([1 for p in g.pprs if p.localization() == 'C']) / total
		loc_m = 100.0 * len([1 for p in g.pprs if p.localization() == 'M']) / total
		loc_s	= 100.0 * len([1 for p in g.pprs if p.localization() == 'S']) / total
		loc_o = 100.0 * len([1 for p in g.pprs if p.localization() == '_']) / total
		loc_u = 100.0 * len([1 for p in g.pprs if p.localization() == '*']) / total
		data.append((short_name(g.name),loc_c,loc_m,loc_s,loc_o,loc_u,))

	data.sort(key=lambda d: d[1])
	utils.write_data(("genome", "c","m","s","other","unknown",),
			data, "output/ppr_localization.dat")

def save_predotar():
	genomes = ppr.get_genomes()
	data = []
	vals =[u'none', u'mitochondrial', u'plastid', u'er', u'elsewhere',
			u'possibly mitochondrial', u'possibly plastid', u'possibly er', 
			u'possibly elsewhere',]
	for g in genomes:
		pprs = ppr.load_records(g)
		if len(pprs) < 50:
			continue

		row = [0,]*len(vals)

		total = float(len(pprs))

		for p in pprs:
			pred = p.annotations['predotar']
			if pred not in vals:
				raise ValueError("didn't expect {}".format(pred))
			else:
				row[vals.index(pred)] += 1
		row = [float(r)/total for r in row]
		data.append([short_name(g),] + row)

	data.sort(key=lambda d: d[1])
	utils.write_data(["genome", ]+vals,
			data, "output/ppr_predotar.dat")


def get_family(pprs):
	try:
		classify.classify([p.seq_record for p in pprs])
	except RuntimeError, e:
		print e
		raw_input("Press enter to continue")
		raise(e)
		
	return (len([1 for p in pprs if p.family() == 'DYW']),
					len([1 for p in pprs if p.family() == 'E+']),
					len([1 for p in pprs if p.family() == 'E']),
					len([1 for p in pprs if p.family() == 'PLS']),
					len([1 for p in pprs if p.family() == 'P']),
				)

def get_length_family(genomes):
	pprs = [p for g in genomes for p in g.pprs]
	type_P   = [p for p in pprs if p.family() == 'P'  ]
	type_E   = [p for p in pprs if p.family() == 'E'  ]
	type_Ep  = [p for p in pprs if p.family() == 'E+' ]
	type_DYW = [p for p in pprs if p.family() == 'DYW']
	type_PLS = [p for p in pprs if p.family() == 'PLS']

	P_hist   = length_hist(type_P)
	E_hist   = length_hist(type_E)
	Ep_hist  = length_hist(type_Ep)
	DYW_hist = length_hist(type_DYW)
	PLS_hist = length_hist(type_PLS)

	hist = fmt_hist([P_hist, E_hist, Ep_hist, DYW_hist, PLS_hist,])

	utils.write_data(('length','p','e','ep','dyw','pls'), hist, 
			"output/ppr_family_lengths.dat")


def get_family_hdr():
	return ("type_dyw", "type_ep", "type_e","type_pls","type_p",)

def short_name(genome):
	i = genome.find('_')
	return genome[0] + genome[i+1].lower()

def main():
	genomes = ppr.get_genomes()
	hist = []
	numbers = []
	locale = []
	family = []
	for g in genomes:
		print g
		pprs = list(ppr.load_pprs(g))
		numbers.append((short_name(g),) + get_numbers(pprs))
		hist.append(length_hist(pprs))
		locale.append((short_name(g),) + get_localization(pprs))
		family.append((short_name(g),) + get_family(pprs))

	numbers.sort(key=lambda n:n[1])
	family.sort(key=lambda f:sum(f[1:]))
	hist = fmt_hist(hist)

	utils.write_data(('genome',)+get_numbers_hdr(), numbers, "output/ppr_numbers.dat")
	utils.write_data(['length',]+[short_name(g) for g in genomes], 
			hist, "output/ppr_lengths.dat")
	utils.write_data(('genome',)+get_localization_hdr(), locale,
		"output/ppr_localization.dat")
	utils.write_data(('genome',)+get_family_hdr(), family, "output/ppr_families.dat")

if __name__ == "__main__":
	main()
