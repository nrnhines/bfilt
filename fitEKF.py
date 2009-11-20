import math
from myscipy import linalg
import numpy
import fitglobals

def ekf(data, model):
        print 'data', data
	debug = fitglobals.debug
	m0 = model.Initial
	P0 = model.P.InitialCov
	k = 0
	smll = 0
	spll = 0
	if debug:
		f = open('data.txt','w')
		f.write('# f1 f2 PLL\n')

	# Main Filtering Loop
	while k < len(data):
		# Evaluate derivatives for prediction
		Times = model.FitEvents[k][0]
		ObsNum = model.FitEvents[k][1]
                print 'model.Obs', model, model.Obs
                print 'ObsNum', ObsNum
		#print "Times m0", Times, m0
		# DEBUGGING CODE
		mb = m0
		for i in range(1,len(Times)):
			 mb = model.Sys.flow([Times[i-1],Times[i]],mb)
			 # if debug:
			 #  print 'flow', mb
		# Operational Code
		# mb = model.Sys.flow(Times, m0)
		# if debug:
		# print 'final flow', mb, '@ Times =', Times
		Am = model.Sys.Dstate(Times, m0)
		#print 'Dstate', Am
		Wm = model.Sys.Dnoise(Times, m0)
		#print "Dnoise", Wm
		
		# Prediction
		Pb = Wm*Wm.T + Am*P0*Am.T

		# Evaluate derivatives for update
		hh = model.Obs.mean(Times, mb, ObsNum)
		H = model.Obs.Dstate(Times, mb, ObsNum)
		V = model.Obs.Dnoise(Times, mb, ObsNum)
		
		# Update
		# print 'data[k]', data[k]
		# print 'ObsNum', ObsNum
                print Times[-1], hh
		if debug:
		  print 'hh', hh
		# print 'data, hh', data[k] hh		
		v = data[k] - hh
		if debug:
		  print 'v', v
		S = H*Pb*H.T + V*V.T
		K = Pb*H.T*S.I
		P0 = Pb - K*S*K.T
                print 'm0', m0
                print 'mb', mb
                print 'K', K
                print 'v', v
		m0 = mb + K*v
		
		# Likelihood
		f0 = len(v)*math.log(2*math.pi)
		f1 =  math.log(linalg.det(S))
		f2 = (v.T*S.I*v).tolist()[0][0]
		smll += (f0+f1+f2)
		k += 1
		if debug:
			print 'S', S
			print 'factors:', f0, f1, f2
			print 'len v:', len(v)
			f.write('%s %s %s\n' % (f1,f2,f1+f2))
		# mll = math.log(linalg.det(S)) + v.T*S.I*v
		# pll = -(1/2)*(f0+f1+f2)
		# spll += pll
		# print 'diff', spll + (1/2)*smll
		# Keeps running sum of Minus-Log-Likelihood
		# But returns Positive Log Likeliood (for now)
	if debug:
		f.close()
	return -smll/2.0
