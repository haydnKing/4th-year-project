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

def get_localization(pprs):
	loc_c = len([p for p in pprs if p.localization() == 'C'])
	loc_m = len([p for p in pprs if p.localization() == 'M'])
	loc_s	= len([p for p in pprs if p.localization() == 'S'])
	loc_other = len([p for p in pprs if p.localization() == '_'])
	loc_unknown = len([p for p in pprs if p.localization() == '*'])
	return (loc_c,loc_m,loc_s,loc_other,loc_unknown,)

def get_localization_hdr():
	return ("loc_c","loc_m","loc_s","loc_other","loc_uk",)

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

def get_family_hdr():
	return ("type_dyw", "type_ep", "type_e","type_pls","type_p",)

def main():
	genomes = ppr.get_processed_genomes()
	hist = []
	numbers = []
	locale = []
	family = []
	for g in genomes:
		print g
		pprs = list(ppr.load_pprs(g))
		numbers.append((g,) + get_numbers(pprs))
		hist.append(length_hist(pprs))
		locale.append((g,) + get_localization(pprs))
		family.append((g,) + get_family(pprs))

	hist = fmt_hist(hist)

	utils.write_data(('genome',)+get_numbers_hdr(), numbers, "output/ppr_numbers.dat")
	utils.write_data(['length',]+genomes, hist, "output/ppr_lengths.dat")
	utils.write_data(('genome',)+get_localization_hdr(), locale,
		"output/ppr_localization.dat")
	utils.write_data(('genome',)+get_family_hdr(), family, "output/ppr_families.dat")

if __name__ == "__main__":
	main()
