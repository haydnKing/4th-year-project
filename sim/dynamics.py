"""Simulate a reaction network involving PPR binding"""
import numpy as np
from scipy.integrate import odeint
from math import log

class Network:
	"""Represent a PPR interaction network"""
	#constants
	vol = 1.0 #litres

	tr_high = 30 / vol #[transcript]/minute
	tr_low  = 0.03 / vol #[transcript]/minute

	tl = 0.693 / vol #[prot]/(min*transcript)

	rna_halflife = 2
	p_halflife = 10

	rna_decay = log(2) / rna_halflife # min^-1
	p_decay   = log(2) / p_halflife # min^-1

	complex_decay_slow = 0.0693 # min^-1
	complex_decay_fast = 0.9 # min^-1

	kD_ppr = 7.6e-9 # Moles / litre
	k_fwd = 2.0 * tr_high
	k_rev = kD_ppr * k_fwd

	def __init__(self, N, interactions, translation):
		"""Build the interaction network
			Arguments:
				N: number of proteins
				interactions: NxN array -
					I[i][j] = 0  if P_i does not bind with RNA_j
					I[i][j] = 1  if (P_i-RNA_j) complex has high translation, slow decay
					I[i][j] =-1  if (P_i-RNA_j) complex has no translation, fast decay
				translation: N array
					t[i] = 0 if RNA_i raw mRNA is not translated,
					t[i] = 1 if RNA_i raw mRNA is translated
		"""
		self.N = N
		self.I = interactions
		self.T = translation

	def simulate(self, inputs, time=60.0, start=None):
		"""Simulate the network dynamics
			inputs: a list of transcriptions (high/low) for each protein
		"""
		N = self.N
		
		def rate(x, t):
			RNA = x[0:self.N]
			P_RNA = []
			for i in range(N):
				P_RNA.append(x[(i+1)*N: (i+2)*N])
			P = x[(N+1)*N:(N+2)*N]

			dRNA = [0.0,]*N
			dP_RNA = [[0.0,]*N,]*N
			dP = [0.0,]*N

			for j in range(N):
				#Transcription
				if inputs[j] > 0:
					dRNA[j] += self.tr_high
				else:
					dRNA[j] += self.tr_low

				#translation
				if self.T[j] > 0:
					dP[j] += RNA[j] * self.tl

				#degredation
				dRNA[j] -= RNA[j] * self.rna_decay
				dP[j] -= P[j] * self.p_decay

				#P_RNA complex formation and decay
				for k in range(N):
					#if P_k binds to RNA_j
					if self.I[k][j] != 0:
						rate = RNA[j]*P[k]*self.k_fwd - P_RNA[k][j] * self.k_rev
						dRNA[j] += -rate 						
						dP[k] +=	-rate
						dP_RNA[k][j] += rate
						#if the interaction allows translation
						if self.I[k][j] > 0:
							dP[j] += P_RNA[k][j] * self.tl
							dP_RNA[k][j] -= P_RNA[k][j] * self.complex_decay_slow
						elif self.I[k][j] < 0:
							dP_RNA[k][j] -= P_RNA[k][j] * self.complex_decay_fast
			
			ret = dRNA
			for i in range(N):
				ret += dP_RNA[i]
			ret += dP
			return ret

		if not start:
			start = [0.0,]*(self.N+2)*self.N
		else:
			p_start = start
			start = [0.0,]*(self.N+1)*self.N + p_start
		t = np.linspace(0, time, 1000*time/60)
		soln = odeint(rate, start, t)

		RNA = soln[:, 0:N]
		P = soln[:, (N+1)*N:(N+2)*N]

		return [t, {'rna': RNA, 'protein': P,},]


