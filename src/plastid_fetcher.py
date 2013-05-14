#!/usr/bin/python

"""Download plastids from entrez"""

import urllib2
from bs4 import BeautifulSoup
from Bio import Entrez
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

def main():
	print "Fetching manifest..."
	lst = get_list()

	for i,(acc, eid) in enumerate(lst):
		stdout.write("\rDownloading ({}/{}): {}         ".format(i,len(lst), acc))
		stdout.flush()
		download(acc, eid)

if __name__ == '__main__':
	main()
