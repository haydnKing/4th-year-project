"""Simulate various logic set-ups"""

import dynamics
from matplotlib import pyplot as plt

def simulate(N, interactions, translation, inputs, outputs, time=60.0):
	"""Simulate and retuen outputs
		N, interactions & translation defined in dynamics
		inputs: list of integers saying which should be treated as inputs -- will
						be simulated as both low and high
		outputs: which protein levels to return"""

	net = dynamics.Network(N, interactions, translation)

	simulations = 2**len(inputs)
	ip = len(inputs)
	sim_inputs = [1,]*N
	fmt = '{{:0{}b}}'.format(ip)

	results = {}

	for sim in range(simulations):
		b = fmt.format(sim)
		for n,i in enumerate(inputs):
			sim_inputs[i] = int(b[n])

		t,out = net.simulate(sim_inputs, time)
		P = out['protein']
		results[b] = [P[:,o] for o in outputs]

	return [t,results]

def NOT():
	out = simulate(2, [[0,-1,],[0,0]], [1,1,], [0,], [1,])
	write_time_sequence(out, "logic_NOT.dat")

def OR():
	out = simulate(3, [[0,0,1,],[0,0,1],[0,0,0,]], [1,1,0,], [0,1,], [2,])
	write_time_sequence(out, "logic_OR.dat")

def NOR():
	out = simulate(3, [[0,0,-1,],[0,0,-1],[0,0,0,]], [1,1,1,], [0,1,], [2,])
	write_time_sequence(out, "logic_NOR.dat")

def NAND():
	out = simulate(5, 
			[[0,0,-1,0,0],[0,0,0,-1,0],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,0],], 
			[1,1,1,1,0], [0,1,], [4,])
	write_time_sequence(out, "logic_NAND.dat")

def Tabor():
	out = simulate(5, 
			[[0,0,0,1,0],[0,0,-1,0,0],[0,0,0,1,0],[0,0,0,0,-1],[0,0,0,0,0],], 
			[1,1,1,0,1], [0,1,], [4,])
	write_time_sequence(out, "logic_Tabor.dat")

def oscillate():
	out = simulate(2,
			[[0,-1,],[1,0,],],
			[0,1,], [], [1,1], 120.0)
	plot(out)

def all():
	NOT()
	OR()
	NOR()
	NAND()
	Tabor()

def plot(out, model='Unnamed'):
	plt.figure()
	(t,results) = out
	for key,value in results.iteritems():
		for v in value:
			plt.plot(t, v, label='{}'.format(key))
	
	plt.xlabel('Time (min)')
	plt.ylabel('Concentration')
	plt.title('Evolution of Protein concentration for {}'.format(model))
	plt.legend(loc=0)
	plt.show()

def write_time_sequence(out, filename):
	t,results = out
	keys = results.keys()
	data = []
	for i,time in enumerate(t):
		data.append((float(time),) + tuple((float(results[key][0][i]) for key in keys)))
	write_data(['time',]+keys, data, filename)

def write_data(header, data, filename):
	"""
		Write data to file in PGF plots friendly formatting
		header = column headers
		data = list of tuples (rows)
		filename = file to write to
	"""
	#Data and header must be present
	if not data:
		raise ValueError("No data given")
	if not header:
		raise ValueError("Header empty")

	#each data item and the header must have the same number of items
	fields = len(header)
	for i,d in enumerate(data, 1):
		if len(d) != fields:
			raise ValueError("Data row {} has {} too {} fields".format(
				i, abs(len(d)-fields), "many" if len(d)>fields else "few"))
	
	#each column must be entirely the same type
	types = tuple((type(d) for d in data[0]))
	for i,d in enumerate(data,1):
		for j,item in enumerate(d):
			if not isinstance(item, types[j]):
				raise TypeError("\'{}\' in row {} is of type {} not {}".format(
					header[j], i, type(item), types[j]))

	formats = {
			str: "{{:<{}s}}",
			int: "{{:>{}d}}",
			float: "{{:>{}.4f}}",
	}
	hformats = {
			str: "{{:<{}s}}",
			int: "{{:>{}s}}",
			float: "{{:>{}s}}",
	}

	#Calculate the maximum width of each
	widths = []
	for i in range(len(header)):
		w = len(formats[str].format(header[i]))
		widths.append(max(w, 
			*[len(formats[types[i]].format(1).format(d[i])) for d in data]))

	try:
		header_fmt = (" "
				.join([hformats[t].format(w) for t,w in zip(types,widths)]) + '\n')
		row_fmt = (" "
				.join([formats[t].format(w) for t,w in zip(types,widths)]) + '\n')
	except KeyError, e:
		raise ValueError("Unsupported type {}".format(e.message))

	f = open(filename, 'w')
	f.write(header_fmt.format(*header))
	for d in data:
		f.write(row_fmt.format(*d))

	f.close()

	return len(data)
