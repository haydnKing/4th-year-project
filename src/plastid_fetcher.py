#!/usr/bin/python

"""Download plastids from entrez"""

import urllib2, os
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
from sys import stdout

LIST_URL = ("http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?"+
														"opt=plastid&amp;taxid=2759&amp;sort=Genome")

Entrez.email = "hjk38@cam.ac.uk"

def get_list():
	"""Return a list of tuples for plastids
		[(accession, entrez_id), ...]
	"""
	f = urllib2.urlopen(LIST_URL)
	soup = BeautifulSoup(f)

	div = soup.find("div", id="Tbl1")
	tbl = div.find("table", "tblTxt")

	ret = []
	for row in tbl.find_all("tr"):
		if "row1" in row.get("class", []):
			continue
		td = row.find_all("td")
		if len(td) == 7:
			link = td[1].find("a")
			try:
				l = link.get("href")
				ret.append( (link.string, l[l.rfind("/")+1:]))
			except AttributeError:
				pass

	return ret

def download(accession, eid):
	f = open("Plastids/{}.gb".format(accession), "w")
	try:
		hnd = Entrez.efetch(db="nucleotide", id=eid, rettype="gb", retmode="text")
		f.write(hnd.read())
	except Exception, e:
		print "Download Failed for \'{}\': {}".format(eid, e)
	hnd.close()
	f.close()

def rename_all():
	for f in os.listdir("Plastids/"):
		try:
			rec = SeqIO.read("Plastids/{}".format(f), 'gb')
			name = rec.description
			cut = max([name.lower().find(s) for s in ['chloroplast', 'apicoplast',
				'plastid', 'chromatophore', 'cyanelle',]])
			if cut < 0:
				print "Failed to find anything for \'{}\'".format(name)
			name = name[0:cut].strip()
			if not name:
				raise ValueError(
					"ERROR: failed to parse \'{}\' from \'{}\'".format(rec.description, f))
			os.remove("Plastids/{}".format(f))
			SeqIO.write(rec, "Plastids/{}.gb".format(name.replace(' ','_')), 'gb')
		except ValueError:
			print "Error in file \'{}\'".format(f)

def main():
	print "Fetching manifest..."
	lst = get_list()

	for i,(acc, eid) in enumerate(lst):
		stdout.write("\rDownloading ({}/{}): {}         ".format(i,len(lst), acc))
		stdout.flush()
		download(acc, eid)

if __name__ == '__main__':
	main()
