#!/usr/bin/python

import subprocess, re, os
from sys import stdout

ROOT_URL = "ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/plants/"
TAIL = "Primary_Assembly/assembled_chromosomes/FASTA/"

chr_re = "^chr(?P<number>\d+).fa.gz$"

def curl(url, outfile=None):
	"""Call curl"""

	if not outfile:
		args = ['curl', '-sl', url,]
	else:
		args = ['curl', '-#', '--output', outfile, url]

	p = subprocess.Popen(args, stdout=subprocess.PIPE)
	stdout, stderr = p.communicate()

	if p.returncode < 0:
		raise ValueError("curl error: {}".format(stderr))

	if not outfile:
		return stdout.split()

def ls(url, exp=None):
	files = curl(url)
	if exp:
		return [f for f in files if re.match(exp, f)]
	return files

def fetch(plant, version, chromasomes):
	head = "{}{}/{}/{}".format(ROOT_URL,plant,version,TAIL)
	outfiles = []
	for i,c in enumerate(chromasomes, 1):
		url = head + c
		outfile = "Genomes/{}".format(c)
		print "\tChromasome {}/{}".format(i, len(chromasomes))
		curl(url, outfile=outfile)
		outfiles.append(outfile)
	
	print "Inflate and merge..."
	outfiles = [inflate(f) for f in outfiles]
	merge(outfiles, "Genomes/{}.fasta".format(plant))

def inflate(f):
	subprocess.call(["gzip", "-d", f])
	m = re.match("^(.*)\.gz", f)
	if m:
		return m.group(1)
	else:
		return None

def merge(in_files, out_file):
	out_file = open(out_file, "w")
	for in_file in in_files:
		i_f = open(in_file, "r")
		data = i_f.read()
		i_f.close()
		out_file.write(data)
		os.remove(in_file)
	out_file.close()

def main():
	stdout.write("\rFetching metadata...")
	stdout.flush()

	plants = ls(ROOT_URL)

	genomes = {}

	for p in plants:
		stdout.write("\rFetching metadata for {}	                       			"
				.format(p))
		stdout.flush()

		plant_url = "{}{}/".format(ROOT_URL, p)
		versions = ls(plant_url)
		genomes[p] = {}
		for v in versions:
			stdout.write("\rFetching metadata for {}/{}	             			      "
				.format(p,v))
			stdout.flush()
			
			version_url = "{}{}/{}".format(plant_url, v, TAIL)
			chromasomes = ls(version_url, exp=chr_re)
			genomes[p][v] = chromasomes

	stdout.write("\r" + (" "*80) + "\r")	
	stdout.flush()

	fetch_list = []

	for plant,versions in genomes.iteritems():

		v_names = [name for name,c in versions.iteritems() if len(c)]
		if not v_names:
			continue
		if len(v_names) == 1:
			selection = v_names[0]
		else:
			print plant
			selection = get_version(versions)
			if not selection:
				continue

		versions[selection].sort(key=lambda c: 
				int(re.match(chr_re,c).group("number")))
		fetch_list.append((plant, selection, versions[selection],))

	print "Fetching Genomes"
	for i,f in enumerate(fetch_list):
		print "({}/{}) {}".format(i, len(fetch_list), f[0])
		fetch(*f)

def get_version(versions):
	v_list = [v for v in versions.keys() if versions[v]]
	print "\t0) None"
	for i,v in enumerate(v_list, 1):
		print "\t{}) {} ({} chr)".format(i, v, len(versions[v]))
	selection = -1
	while selection not in range(len(v_list)+1):
		try:
			selection = int(raw_input("Select a version by number: "))
		except ValueError:
			selection = -1
	
	if selection == 0:
		return None		
	return v_list[selection -1]



if __name__ == '__main__':
	main()
