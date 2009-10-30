import math
from myscipy import linalg
import numpy

def ekf(data, model):
	m0 = model.Initial
	P0 = model.P.InitialCov
	k = 0
	smll = 0
	
	# Main Filtering Loop
	while k < len(data):
		# Evaluate derivatives for prediction
		Times = model.FitEvents[k][0]
		ObsNum = model.FitEvents[k][1]
		#print "Times m0", Times, m0
		# DEBUGGING CODE
		mb = m0
		for i in range(1,len(Times)):
			 mb = model.Sys.flow([Times[i-1],Times[i]],mb)
			 print 'flow', mb
		# Operational Code
		# mb = model.Sys.flow(Times, m0)
		print 'final flow', mb, '@ Times =', Times
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
		print 'hh', hh
		# print 'data, hh', data[k] hh		
		v = data[k] - hh
		print 'v', v
		S = H*Pb*H.T + V*V.T
		K = Pb*H.T*S.I
		P0 = Pb - K*S*K.T
		m0 = mb + K*v
		
		# Likelihood
		mll = math.log(linalg.det(S)) + v.T*S.I*v
		smll += mll
		k += 1
		# Keeps running sum of Minus-Log-Likelihood
		# But returns Positive Log Likeliood (for now)
	return -smll
