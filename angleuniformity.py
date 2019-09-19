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
	theta = np.arctan2(y,x)
	phi = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return phi, theta, r

thetas=[0.5625
,1.6875
,2.8125
,3.9375
,5.0625
,6.1875
,7.3125
,8.4375
,9.5625
,10.6875
,11.8125
,12.9375
,14.0625
,15.1875
,16.3125
,17.4375
,18.5625
,19.6875
,20.8125
,21.9375
,23.0625
,24.1875
,25.3125
,26.4375
,27.5625
,28.6875
,29.8125
,30.9375
,32.0625
,33.1875
,34.3125
,35.4375
,36.5625
,37.6875
,38.8125
,39.9375
,41.0625
,42.1875
,43.3125
,44.4375
,45.5625
,46.6875
,47.8125
,48.9375
,50.0625
,51.1875
,52.3125
,53.4375
,54.5625
,55.6875
,56.8125
,57.9375
,59.0625
,60.1875
,61.3125
,62.4375
,63.5625
,64.6875
,65.8125
,66.9375
,68.0625
,69.1875
,70.3125
,71.4375
,72.5625
,73.6875
,74.8125
,75.9375
,77.0625
,78.1875
,79.3125
,81.5625
,82.6875]

outputfile = TFile("ElectronAngleUniformity.root", "RECREATE")
for theta in thetas:

	file = "/home/software/Calo/results/Thetauniformity/Theta_"+str(theta)+".txt"	

	if theta == thetas[0]:
		truethetas = array('d')
		angrestheta = array('d')
		thetas = array('d')

	EvtID = np.array([line.split('\t')[0] for line in open(file,"r").readlines()][1:],'i')
	NumFiber = np.array([x.split('\t')[1] for x in open(file,"r").readlines()][1:],'d')
	Edep = np.array([x.split('\t')[2] for x in open(file,"r").readlines()][1:],'d')
	X = np.array([x.split('\t')[3] for x in open(file,"r").readlines()][1:],'d')
	Y = np.array([x.split('\t')[4] for x in open(file,"r").readlines()][1:],'d')
	Z = np.array([x.split('\t')[5] for x in open(file,"r").readlines()][1:],'d')
	Flag = np.array([x.split('\t')[6] for x in open(file,"r").readlines()][1:],'d')
	Slice = np.array([x.split('\t')[7] for x in open(file,"r").readlines()][1:],'d')
	Tower = np.array([x.split('\t')[8] for x in open(file,"r").readlines()][1:],'d')	

	thetarad = theta*math.pi/180.

	ThetaHist = TH1F("Theta"+str(theta), "Theta", 500, thetarad-0.01, thetarad+0.01)
	PhiHist = TH1F("Phi"+str(theta), "Phi", 500, -0.01, 0.01)	

	ThetaHistC = TH1F("Theta_C"+str(theta), "Theta", 500, thetarad-0.01, thetarad+0.01)
	PhiHistC = TH1F("Phi_C"+str(theta), "Phi", 500, -0.01, 0.01)	

	ThetaHistCS = TH1F("Theta_CS"+str(theta), "Theta", 500, thetarad-0.01, thetarad+0.01)
	PhiHistCS = TH1F("Phi_CS"+str(theta), "Phi", 500, -0.01, 0.01)	

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

		ScinPlot = TH2F("Scin_"+str(theta), "Scin_"+str(theta), int(100), float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		CherPlot = TH2F("Cher_"+str(theta), "Cher_"+str(theta), int(100),  float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		
		n = len(theta_S)
		array_theta_S = np.array(theta_S, 'd')
		array_phi_S = np.array(phi_S, 'd')
		array_E_s = np.array(E_s, 'd')
		ScinGraph = TGraph2D(n, array_theta_S, array_phi_S, array_E_s)
		ScinGraph.SetName("ScinGraph_"+str(theta))

		n = len(theta_C)
		array_theta_C = np.array(theta_C, 'd')
		array_phi_C = np.array(phi_C, 'd')
		array_E_c = np.array(E_c, 'd')
		CherGraph = TGraph2D(n, array_theta_C, array_phi_C, array_E_c)
		CherGraph.SetName("CherGraph_"+str(theta))

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

	truethetas.append(theta)
	angrestheta.append(ThetaHistCS.GetFunction("gaus").GetParameter(2)*1000)
	thetas.append(ThetaHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)

	print "true theta "+str(theta)+" measured theta "+str(ThetaHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)+" sigma theta "+str(ThetaHistCS.GetFunction("gaus").GetParameter(2)*1000)

gStyle.SetOptStat(111)
theta_linearity = TGraph(len(truethetas), truethetas, thetas)
theta_linearity.SetName("theta_linearity")
theta_linearity.Write()
sigma_linearity = TGraph(len(truethetas), truethetas, angrestheta)
sigma_linearity.SetName("sigma_linearity")
sigma_linearity.Write()