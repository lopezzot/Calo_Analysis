import os
import numpy as np
import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F
import circle

def cart2sph(x,y,z):
    theta = np.arctan2(y,x)
    phi = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    return theta, phi, r

outputfile = TFile("Displays.root", "RECREATE")

#file = "/Users/lorenzo/Desktop/Git_IDEA_CALO_FIBER-build/B4a/Event0-0.txt"
file = "Event-0-0-40GeV.txt"

EvtID = np.array([line.split('\t')[0] for line in open(file,"r").readlines()][1:],'i')
NumFiber = np.array([x.split('\t')[1] for x in open(file,"r").readlines()][1:],'d')
Edep = np.array([x.split('\t')[2] for x in open(file,"r").readlines()][1:],'d')
X = np.array([x.split('\t')[3] for x in open(file,"r").readlines()][1:],'d')
Y = np.array([x.split('\t')[4] for x in open(file,"r").readlines()][1:],'d')
Z = np.array([x.split('\t')[5] for x in open(file,"r").readlines()][1:],'d')
Flag = np.array([x.split('\t')[6] for x in open(file,"r").readlines()][1:],'d')
Slice = np.array([x.split('\t')[7] for x in open(file,"r").readlines()][1:],'d')
Tower = np.array([x.split('\t')[8] for x in open(file,"r").readlines()][1:],'d')

ThetaHist = TH1F("Theta", "Theta", 1000, -0.02, 0.02)
PhiHist = TH1F("Phi", "Phi", 1000, -0.02, 0.02)
percentages_array = [0.0,0.0,0.0,0.0,0.0]
percentages_array = np.array(percentages_array)
for i in range (0,int(max(EvtID[1:]))+1):
	# S have Flag == 1 C == 0
	S = np.where((Flag==1.) & (EvtID==i))
	X_S=X[S]
	Y_S=Y[S]
	Z_S=Z[S]
	theta_S=[]
	phi_S=[]
	E_s = Edep[S]
	C = np.where((Flag==0.) & (EvtID==i))
	X_C=X[C]
	Y_C=Y[C]
	Z_C=Z[C]
	theta_C=[]
	phi_C=[]
	E_c=Edep[C]

	for counter, j in enumerate(X_S):
		(Stheta, Sphi, Sr)=cart2sph(float(X_S[counter]),float(Y_S[counter]),float(Z_S[counter]))
		theta_S.append(Stheta)
		phi_S.append(Sphi)	
	
	for counter, j in enumerate(X_C):
		(Ctheta, Cphi, Cr)=cart2sph(float(X_C[counter]),float(Y_C[counter]),float(Z_C[counter]))
		theta_C.append(Ctheta)
		phi_C.append(Cphi)	

	
	ScinPlot = TH2F("Scin_"+str(i), "Scin_"+str(i), int(100), float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
	
	MeanTheta = 0.
	MeanPhi = 0.
	sumtheta = 0.
	sumphi = 0.

	for counter, entry in enumerate(E_s):
		MeanTheta += entry*theta_S[counter]
		sumtheta += entry
		MeanPhi += entry*phi_S[counter]
		sumphi += entry
		ScinPlot.Fill(theta_S[counter], phi_S[counter], entry)

	MeanTheta = MeanTheta/sumtheta
	MeanPhi = MeanPhi/sumphi
	ThetaHist.Fill(MeanTheta)
	PhiHist.Fill(MeanPhi)
	percentages = circle.computerradi(E_s, theta_S, phi_S)
	
	percentages_array = percentages_array + np.array(percentages) 
	print str(i)+" "+str(percentages)
	if i<10:
		ScinPlot.Write()
	del ScinPlot

	#plt.plot(theta_S, phi_S,'.b',ms=1)
	#plt.show()
	
	'''
	for h in np.arange(0,len(C[0][:])):
		(Ctheta, Cphi, Cr)=cart2sph(float(X_C[h]),float(Y_C[h]),float(Z_C[h]))
		theta_C.append(Ctheta)
		phi_C.append(Cphi)	
	''' 
	#plt.plot(theta_C, phi_C,'.r',ms=1)
	#plt.show()

percentages_array = percentages_array/3000
print percentages_array
ThetaHist.Write()
PhiHist.Write()
