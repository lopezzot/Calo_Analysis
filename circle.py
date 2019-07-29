

def computerradi(E_s, theta_S, phi_S):

	#compute bariocenter
	MeanTheta = 0.
	MeanPhi = 0.
	sumtheta = 0.
	sumphi = 0.

	for counter, entry in enumerate(E_s):
		MeanTheta += entry*theta_S[counter]
		sumtheta += entry
		MeanPhi += entry*phi_S[counter]
		sumphi += entry

	MeanTheta = MeanTheta/sumtheta
	MeanPhi = MeanPhi/sumphi
	#end bariocenter

	radii = [0.005, 0.01, 0.05, 0.1, 0.15]
	percentages = []
	for radius in radii:
		percentages.append(computerradius(E_s, theta_S, phi_S, MeanTheta, MeanPhi, radius))
		
	return percentages


def computerradius(E_s, Theta_s, Phi_s, Theta_b, Phi_b, radius):
	'''function to perform energy inside circle'''
	signalcontained = 0.
	for counter, signal in enumerate(E_s):
		if ((Theta_s[counter]-Theta_b)**2+(Phi_s[counter]-Phi_b)**2)**0.5< radius:
			#print ((Theta_s[counter]-Theta_b)**2+(Phi_s[counter]-Phi_b)**2)**0.5
			signalcontained += signal

	return signalcontained/sum(E_s)

