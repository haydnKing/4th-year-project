"""Call targetp and annotate the SeqRecord(s) with the result"""

import subprocess, cStringIO as StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

location_string = {
						'C': 'Chloroplast',
						'M': 'Mitochondrion',
						'S': 'Secretory',
						'_': 'Other',
						'*': 'Unknown',
						}


def targetp(recs, annotation='target'):
	"""annotate the record(s) with their predicted location
		rec can be a single SeqRecord or a list of them
	"""
	if isinstance(recs, SeqRecord):
		recs = [recs,]
		
	#replace id and desc to be distinct
	data = []
	for i,r in enumerate(recs):
		data.append((r.id, r.description, r,))
		r.id = str(i)
		r.description = ''
		print "r.id = {}".format(r.id)

	s = StringIO.StringIO()
	SeqIO.write(recs, s, 'fasta')

	p = subprocess.Popen(['targetp', '-P',], stdin=subprocess.PIPE, 
			stdout=subprocess.PIPE)
	(stdout, stderr) = p.communicate(input=s.getvalue())

	state = False
	print stdout
	for line in stdout.splitlines():
		if not state:
			if line.find('--------') >= 0:
				state = True
				continue
		else:
			if line.find('--------') >= 0:
				state = False
				continue
			
			l = line.split()
			try:
				data[int(l[0])][2].annotations[annotation] = l[6]
			except ValueError:
				continue

	
	#put recs back as it was
	for d in data:
		d[2].id = d[0]
		d[2].description = d[1]
	
