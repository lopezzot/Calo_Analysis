import os
import numpy as np
#import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F, TGraph2D, TLine, gStyle, TGraph
import circle
import createangularesolutionplot
from array import array
import math

def cart2sph(x,y,z):
	if y > 0.:
		phi = np.arccos(x/np.sqrt(x**2+y**2))
	else:
		phi = 2*math.pi-np.arccos(x/np.sqrt(x**2+y**2))
	theta = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return phi, theta, r

phis=[5.0,
15.0,
25.0,
35.0,
45.0,
55.0,
65.0,
75.0,
85.0,
95.0,
105.0,
115.0,
125.0,
135.0,
145.0,
155.0,
165.0,
175.0]
'''
185.0,
195.0,
205.0,
215.0,
225.0,
235.0,
245.0,
255.0,
265.0,
275.0,
285.0,
295.0,
305.0,
315.0,
325.0,
335.0,
345.0,
355.0]
'''
phis = [25.0]

outputfile = TFile("ElectronAngleUniformity_phi.root", "RECREATE")
for phi in phis:

	file = "/home/software/Calo/results/philinearity/Phi_"+str(phi)+".txt"	
	file = "/home/lorenzo/Desktop/Calo/results/philinearity/Phi_"+str(phi)+".txt"	

	if phi == phis[0]:
		truephis = array('d')
		angresphi = array('d')
		phis = array('d')

	EvtID = np.array([line.split('\t')[0] for line in open(file,"r").readlines()][1:],'i')
	NumFiber = np.array([x.split('\t')[1] for x in open(file,"r").readlines()][1:],'d')
	Edep = np.array([x.split('\t')[2] for x in open(file,"r").readlines()][1:],'d')
	X = np.array([x.split('\t')[3] for x in open(file,"r").readlines()][1:],'d')
	Y = np.array([x.split('\t')[4] for x in open(file,"r").readlines()][1:],'d')
	Z = np.array([x.split('\t')[5] for x in open(file,"r").readlines()][1:],'d')
	Flag = np.array([x.split('\t')[6] for x in open(file,"r").readlines()][1:],'d')
	Slice = np.array([x.split('\t')[7] for x in open(file,"r").readlines()][1:],'d')
	Tower = np.array([x.split('\t')[8] for x in open(file,"r").readlines()][1:],'d')	

	phirad = phi*math.pi/180.

	ThetaHist = TH1F("Theta"+str(phi), "phi", 500, -0.01, 0.01)
	PhiHist = TH1F("Phi"+str(phi), "Phi", 500, phirad-0.01, phirad+0.01)	

	ThetaHistC = TH1F("Theta_C"+str(phi), "phi", 500, -0.01, 0.01)
	PhiHistC = TH1F("Phi_C"+str(phi), "Phi", 500, phirad-0.01, phirad+0.01)	

	ThetaHistCS = TH1F("Theta_CS"+str(phi), "phi", 500, -0.01, 0.01)
	PhiHistCS = TH1F("Phi_CS"+str(phi), "Phi", 500, phirad-0.01, phirad+0.01)	

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
			(Sphi, Stheta, Sr)=cart2sph(float(X_S[counter]),float(Y_S[counter]),float(Z_S[counter]))
			print float(X_S[counter]),float(Y_S[counter]),float(Z_S[counter]),Sphi
			theta_S.append(Stheta)
			'''
			if phi >= 185.0:
				Sphi = math.pi*2-abs(Sphi)
				phi_S.append(Sphi)
			'''
			phi_S.append(Sphi)	
		
		for counter, j in enumerate(X_C):
			(Cphi, Ctheta, Cr)=cart2sph(float(X_C[counter]),float(Y_C[counter]),float(Z_C[counter]))
			theta_C.append(Ctheta)
			'''
			if phi >= 185.0:
				Cphi = math.pi*2-abs(Cphi)
				phi_C.append(Cphi)
			'''			
			phi_C.append(Cphi)		
		
		ScinPlot = TH2F("Scin_"+str(phi), "Scin_"+str(phi), int(100), float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		CherPlot = TH2F("Cher_"+str(phi), "Cher_"+str(phi), int(100),  float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		
		n = len(theta_S)
		array_theta_S = np.array(theta_S, 'd')
		array_phi_S = np.array(phi_S, 'd')
		array_E_s = np.array(E_s, 'd')
		ScinGraph = TGraph2D(n, array_theta_S, array_phi_S, array_E_s)
		ScinGraph.SetName("ScinGraph_"+str(phi))

		n = len(theta_C)
		array_theta_C = np.array(theta_C, 'd')
		array_phi_C = np.array(phi_C, 'd')
		array_E_c = np.array(E_c, 'd')
		CherGraph = TGraph2D(n, array_theta_C, array_phi_C, array_E_c)
		CherGraph.SetName("CherGraph_"+str(phi))

		MeanTheta = 0.
		MeanPhi = 0.
		sumtheta = 0.
		sumphi = 0.
		MeanThetaC = 0.
		MeanPhiC = 0.
		sumthetaC = 0.
		sumphiC = 0.
		MeanThetaCS = 0.
		MeanPhiCS = 0.
		sumthetaCS = 0.
		sumphiCS = 0.	

		for counter, entry in enumerate(E_s):
			MeanTheta += entry*theta_S[counter]
			sumtheta += entry
			MeanPhi += entry*phi_S[counter]
			sumphi += entry
			ScinPlot.Fill(theta_S[counter], phi_S[counter], entry)	

		for counter, entry in enumerate(E_c):
			MeanThetaC += entry*theta_C[counter]
			sumthetaC += entry
			MeanPhiC += entry*phi_C[counter]
			sumphiC += entry
			CherPlot.Fill(theta_C[counter], phi_C[counter], entry)	

		for counter, entry in enumerate(E_s):
			MeanThetaCS += entry*theta_S[counter]
			sumthetaCS += entry
			MeanPhiCS += entry*phi_S[counter]
			sumphiCS += entry
		for counter, entry in enumerate(E_c):
			MeanThetaCS += entry*theta_C[counter]
			sumthetaCS += entry
			MeanPhiCS += entry*phi_C[counter]
			sumphiCS += entry
			
		MeanTheta = MeanTheta/sumtheta
		MeanPhi = MeanPhi/sumphi
		ThetaHist.Fill(MeanTheta)
		PhiHist.Fill(MeanPhi)
		
		MeanThetaC = MeanThetaC/sumthetaC
		MeanPhiC = MeanPhiC/sumphiC
		ThetaHistC.Fill(MeanThetaC)
		PhiHistC.Fill(MeanPhiC)	

		MeanThetaCS = MeanThetaCS/sumthetaCS
		MeanPhiCS = MeanPhiCS/sumphiCS
		ThetaHistCS.Fill(MeanThetaCS)
		PhiHistCS.Fill(MeanPhiCS)	

		if i<1:
			gStyle.SetPalette(1)
			ScinPlot.Write()
			CherPlot.Write()
			ScinGraph.Write()
			CherGraph.Write()
		
		del ScinPlot
		del CherPlot	
		del ScinGraph
		del CherGraph

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

	ThetaHist.Fit("gaus")
	PhiHist.Fit("gaus")
	ThetaHistC.Fit("gaus")
	PhiHistC.Fit("gaus")
	ThetaHistCS.Fit("gaus")
	PhiHistCS.Fit("gaus")

	ThetaHist.Write()
	PhiHist.Write()
	ThetaHistC.Write()
	PhiHistC.Write()
	ThetaHistCS.Write()
	PhiHistCS.Write()
	truephis.append(phi)
	angresphi.append(PhiHistCS.GetFunction("gaus").GetParameter(2)*1000)
	phis.append(PhiHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)

	print "true phi "+str(phi)+" measured theta "+str(PhiHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)+" sigma phi "+str(PhiHistCS.GetFunction("gaus").GetParameter(2)*1000)

gStyle.SetOptStat(111)
phi_linearity = TGraph(len(truephis), truephis, phis)
phi_linearity.SetName("phis_linearity")
phi_linearity.Write()
sigma_linearity = TGraph(len(truephis), truephis, angresphi)
sigma_linearity.SetName("sigma_linearity")
sigma_linearity.Write()