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
	plot(out, 'NOT')

def OR():
	out = simulate(3, [[0,0,1,],[0,0,1],[0,0,0,]], [1,1,0,], [0,1,], [2,])
	plot(out, 'OR')

def NOR():
	out = simulate(3, [[0,0,-1,],[0,0,-1],[0,0,0,]], [1,1,1,], [0,1,], [2,])
	plot(out, 'NOR')

def NAND():
	out = simulate(5, 
			[[0,0,-1,0,0],[0,0,0,-1,0],[0,0,0,0,1],[0,0,0,0,1],[0,0,0,0,0],], 
			[1,1,1,1,0], [0,1,], [4,])
	plot(out, 'NAND')

def Tabor():
	out = simulate(5, 
			[[0,0,0,1,0],[0,0,-1,0,0],[0,0,0,1,0],[0,0,0,0,-1],[0,0,0,0,0],], 
			[1,1,1,0,1], [0,1,], [4,])
	plot(out, 'Tabor')

def oscillate():
	out = simulate(2,
			[[0,-1,],[1,0,],],
			[0,1,], [], [1,1], 120.0)
	plot(out)

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


