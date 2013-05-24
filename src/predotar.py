import urllib2, urllib
from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

INTS = ['elsewhere', 'er', 'mitochondrial','plastid', ]
url = "http://urgi.versailles.inra.fr/predotar/test.seq"

def parse_response(data):
	soup = BeautifulSoup(data)

	table = soup.find_all("table")[1]

	rows = table.find_all("tr")

	ret = []
	keys = []
	for td in rows[0]:
		keys.append(td.string.lower())

	for row in rows[1:]:
		val = {}
		for key,td in zip(keys, row("td")):
			try:
				value = td.string
				if key in INTS:
					val[key] = float(value.replace(',','.'))
				else:
					val[key] = value
			except ValueError:
				continue
		ret.append(val)
	
	return ret

def get_prediction(pprs):
	f = StringIO()
	
	dpprs = [SeqRecord(p.seq.translate(), "{}".format(i)) for i,p in
			enumerate(pprs)]
	
	SeqIO.write(dpprs, f, "fasta")
	data = urllib.urlencode({'name':'',
														'sequence': f.getvalue(),
														'plastids': 1,})

	output = parse_response(urllib2.urlopen(url=url, data=data).read())

	for o in output:
		i = int(o['sequence'])
		try:
			pprs[i].annotations['predotar'] = o['prediction'].lower()
		except KeyError:
			pprs[i].annotations['predotar'] = u'none'

def predict(p):
	if isinstance(p, SeqRecord):
		get_prediction([p,])
	else:
		get_prediction(p)

	



	
