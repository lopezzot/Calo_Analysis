import os
import numpy as np
#import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F, TGraph2D, TLine, gStyle, TGraph, TVector3
import circle
import createangularesolutionplot
from array import array
import math

def cart2sph(x,y,z):
	if y > 0.:
		phi = np.arccos(x/np.sqrt(x**2+y**2))
	else:
		phi = 2*math.pi-np.arccos(x/np.sqrt(x**2+y**2))
	#phi = np.arctan2(y,x)
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
175.0,
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

outputfile = TFile("ElectronAngleUniformity_phi.root", "RECREATE")
for phi in phis:

	file = "/home/software/Calo/results/Phiuniformity/Phi_"+str(phi)+".txt"	
	#file = "/home/lorenzo/Desktop/Calo/results/philinearity/Phi_"+str(phi)+".txt"	

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

	ThetaHist = TH1F("Theta"+str(phi), "theta", 300, -0.01, 0.01)
	PhiHist = TH1F("Phi"+str(phi), "Phi", 300, phirad-0.01, phirad+0.01)		
	ThetaHistC = TH1F("ThetaC"+str(phi), "theta", 300, -0.01, 0.01)
	PhiHistC = TH1F("PhiC"+str(phi), "Phi", 300, phirad-0.01, phirad+0.01)		
	ThetaHistCS = TH1F("ThetaCS"+str(phi), "theta", 300, -0.01, 0.01)
	PhiHistCS = TH1F("PhiCS"+str(phi), "Phi", 300, phirad-0.01, phirad+0.01)		

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

		scin_xarray = np.array(X_S,'d')
		scin_yarray = np.array(Y_S,'d')
		scin_zarray = np.array(Z_S,'d')
		scin_array = np.array(E_s,'d')

		cher_xarray = np.array(X_C,'d')
		cher_yarray = np.array(Y_C,'d')
		cher_zarray = np.array(Z_C,'d')
		cher_array = np.array(E_c,'d')
		
		MeanX = 0.
		MeanY = 0.
		MeanZ = 0.
		sumscin = 0.

		cher_MeanX = 0.
		cher_MeanY = 0.
		cher_MeanZ = 0.
		sumcher = 0.

		for counter, entry in enumerate(scin_array):
			MeanX += entry*scin_xarray[counter]
			MeanY += entry*scin_yarray[counter]
			MeanZ += entry*scin_zarray[counter]
			sumscin += entry

		for counter, entry in enumerate(cher_array):
			cher_MeanX += entry*cher_xarray[counter]
			cher_MeanY += entry*cher_yarray[counter]
			cher_MeanZ += entry*cher_zarray[counter]
			sumcher += entry
				
		MeanX = MeanX/sumscin
		MeanY = MeanY/sumscin
		MeanZ = MeanZ/sumscin
		
		cher_MeanX = cher_MeanX/sumcher
		cher_MeanY = cher_MeanY/sumcher
		cher_MeanZ = cher_MeanZ/sumcher
		
		M_phi, M_theta, M_r = cart2sph(MeanX, MeanY, MeanZ)
		ThetaHist.Fill(M_theta)
		PhiHist.Fill(M_phi)

		cher_M_phi, cher_M_theta, cher_M_r = cart2sph(cher_MeanX, cher_MeanY, cher_MeanZ)
		ThetaHistC.Fill(cher_M_theta)
		PhiHistC.Fill(cher_M_phi)

		ThetaHistCS.Fill((cher_M_theta+M_theta)/2.)
		PhiHistCS.Fill((cher_M_phi+M_phi)/2.)

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
	phis.append(PhiHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi-phi)

	print "true phi "+str(phi)+" measured phi "+str(PhiHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)+" sigma phi "+str(PhiHistCS.GetFunction("gaus").GetParameter(2)*1000)

gStyle.SetOptStat(111)
phi_linearity = TGraph(len(truephis), truephis, phis)
phi_linearity.SetName("phis_linearity")
phi_linearity.Write()
sigma_linearity = TGraph(len(truephis), truephis, angresphi)
sigma_linearity.SetName("sigma_linearity")
sigma_linearity.Write()