import os
import numpy as np
#import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F, TGraph2D, TLine, gStyle, TGraph
import circle
import createangularesolutionplot

def cart2sph(x,y,z):
	theta = np.arctan2(y,x)
	phi = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return phi, theta, r
energies = [30,40,50,60,70,80,90,100,110,120,130,140,150]

outputfile = TFile("PionAngle.root", "RECREATE")
for e in energies:

	#file = "/Users/lorenzo/Desktop/Git_IDEA_CALO_FIBER-build/B4a/Event0-0.txt"
	#file = "Event-0-0-40GeV.txt"
	file = "/home/software/Calo/results/energy_angular_pion/Pion_"+str(e)+".txt"	
	#file = "/home/software/Calo/results/angular_res/Electron_"+str(e)+".txt"

	if e == energies[0]:
		energies = []
		angrestheta = []
		angresphi = []
		angrestheta_c = []
		angrestheta_s = []
		angresphi_c = []
		angresphi_s = []

	EvtID = np.array([line.split('\t')[0] for line in open(file,"r").readlines()][1:],'i')
	NumFiber = np.array([x.split('\t')[1] for x in open(file,"r").readlines()][1:],'d')
	Edep = np.array([x.split('\t')[2] for x in open(file,"r").readlines()][1:],'d')
	X = np.array([x.split('\t')[3] for x in open(file,"r").readlines()][1:],'d')
	Y = np.array([x.split('\t')[4] for x in open(file,"r").readlines()][1:],'d')
	Z = np.array([x.split('\t')[5] for x in open(file,"r").readlines()][1:],'d')
	Flag = np.array([x.split('\t')[6] for x in open(file,"r").readlines()][1:],'d')
	Slice = np.array([x.split('\t')[7] for x in open(file,"r").readlines()][1:],'d')
	Tower = np.array([x.split('\t')[8] for x in open(file,"r").readlines()][1:],'d')	

	ThetaHist = TH1F("Theta"+str(e), "Theta", 500, -0.03, 0.03)
	PhiHist = TH1F("Phi"+str(e), "Phi", 500, -0.03, 0.03)	

	ThetaHistC = TH1F("Theta_C"+str(e), "Theta", 500, -0.03, 0.03)
	PhiHistC = TH1F("Phi_C"+str(e), "Phi", 500, -0.03, 0.03)	

	ThetaHistCS = TH1F("Theta_CS"+str(e), "Theta", 500, -0.03, 0.03)
	PhiHistCS = TH1F("Phi_CS"+str(e), "Phi", 500, -0.03, 0.03)	

	percentages_array = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
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

		ScinPlot = TH2F("Scin_"+str(e), "Scin_"+str(e), int(100), float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		CherPlot = TH2F("Cher_"+str(e), "Cher_"+str(e), int(100),  float(np.mean(theta_S)-0.05), float(np.mean(theta_S)+0.05), int(100), float(np.mean(phi_S)-0.05),float(np.mean(phi_S)+0.05))
		
		n = len(theta_S)
		array_theta_S = np.array(theta_S, 'd')
		array_phi_S = np.array(phi_S, 'd')
		array_E_s = np.array(E_s, 'd')
		ScinGraph = TGraph2D(n, array_theta_S, array_phi_S, array_E_s)
		ScinGraph.SetName("ScinGraph_"+str(e))

		n = len(theta_C)
		array_theta_C = np.array(theta_C, 'd')
		array_phi_C = np.array(phi_C, 'd')
		array_E_c = np.array(E_c, 'd')
		CherGraph = TGraph2D(n, array_theta_C, array_phi_C, array_E_c)
		CherGraph.SetName("CherGraph_"+str(e))

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
		#percentages = circle.computerradi(E_s, theta_S, phi_S)
		
		MeanThetaC = MeanThetaC/sumthetaC
		MeanPhiC = MeanPhiC/sumphiC
		ThetaHistC.Fill(MeanThetaC)
		PhiHistC.Fill(MeanPhiC)	
		percentages = circle.computerradi(E_c, theta_C, phi_C)

		MeanThetaCS = MeanThetaCS/sumthetaCS
		MeanPhiCS = MeanPhiCS/sumphiCS
		ThetaHistCS.Fill(MeanThetaCS)
		PhiHistCS.Fill(MeanPhiCS)	

		percentages_array = percentages_array + np.array(percentages) 
		print str(i)+" "+str(percentages)
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

	percentages_array = percentages_array/3000
	print percentages_array
	radii = np.array([0.001, 0.003, 0.005, 0.01, 0.05, 0.1, 0.15], 'd')
	p = TGraph(len(radii), radii, percentages_array)
	p.Write()

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

	energies.append(float(e))
	angrestheta.append(ThetaHistCS.GetFunction("gaus").GetParameter(2)*1000)
	angresphi.append(PhiHistCS.GetFunction("gaus").GetParameter(2)*1000)
	angrestheta_s.append(ThetaHist.GetFunction("gaus").GetParameter(2)*1000)
	angresphi_s.append(PhiHist.GetFunction("gaus").GetParameter(2)*1000)
	angrestheta_c.append(ThetaHistC.GetFunction("gaus").GetParameter(2)*1000)
	angresphi_c.append(PhiHistC.GetFunction("gaus").GetParameter(2)*1000)
		
	print angrestheta, angresphi, energies
gStyle.SetOptStat(111)
ThetaGraph, PhiGraph = createangularesolutionplot.angresplot(energies, angrestheta, angresphi)
ThetaGraph_s, PhiGraph_s = createangularesolutionplot.angresplot(energies, angrestheta_s, angresphi_s)
ThetaGraph_c, PhiGraph_c = createangularesolutionplot.angresplot(energies, angrestheta_c, angresphi_c)
ThetaGraph.SetName("ThetaGraph")
ThetaGraph.Write()
PhiGraph.SetName("PhiGraph")
PhiGraph.Write()
ThetaGraph_s.SetName("ThetaGraph_s")
ThetaGraph_s.Write()
PhiGraph_s.SetName("PhiGraph_s")
PhiGraph_s.Write()
ThetaGraph_c.SetName("ThetaGraph_c")
ThetaGraph_c.Write()
PhiGraph_c.SetName("PhiGraph_c")
PhiGraph_c.Write()