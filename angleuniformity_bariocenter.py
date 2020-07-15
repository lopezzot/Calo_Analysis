import os
import numpy as np
#import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F, TGraph2D, TLine, gStyle, TGraph, TGraphErrors
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
	#file = "/home/lorenzo/Desktop/Calo/results/thetalinearity/Theta_"+str(theta)+".txt"	
	print file
	if theta == thetas[0]:
		truetheta = array('d')
		angrestheta = array('d')
		angresthetaerror = array('d')
		thetas = array('d')
		truethetaerror = array('d')

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

	ThetaHist = TH1F("Theta"+str(theta), "theta", 300, thetarad-0.01, thetarad+0.01)
	PhiHist = TH1F("Phi"+str(theta), "Phi", 300,  -0.01, 0.01)		
	ThetaHistC = TH1F("ThetaC"+str(theta), "theta", 300, thetarad-0.01, thetarad+0.01)
	PhiHistC = TH1F("PhiC"+str(theta), "Phi", 300,  -0.01, 0.01)		
	ThetaHistCS = TH1F("ThetaCS"+str(theta), "theta", 300, thetarad-0.01, thetarad+0.01)
	PhiHistCS = TH1F("PhiCS"+str(theta), "Phi", 300, -0.01, 0.01)		

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

	truetheta.append(90.+theta)
	truethetaerror.append(0.)
	angrestheta.append(ThetaHistCS.GetFunction("gaus").GetParameter(2)*1000)
	angresthetaerror.append(ThetaHistCS.GetFunction("gaus").GetParError(2)*1000)
	thetas.append(90.+ThetaHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi-(90.+theta))

	print "true theta "+str(theta)+" measured theta "+str(ThetaHistCS.GetFunction("gaus").GetParameter(1)*180./math.pi)+" sigma theta "+str(ThetaHistCS.GetFunction("gaus").GetParameter(2)*1000)

gStyle.SetOptStat(111)
phi_linearity = TGraph(len(truetheta), truetheta, thetas)
phi_linearity.SetName("thetas_linearity")
phi_linearity.GetXaxis().SetTitle("#theta (deg)")
phi_linearity.GetYaxis().SetTitle("#theta measured - theta true (deg)")
phi_linearity.Write()
sigma_linearity = TGraphErrors(len(truetheta), truetheta, angrestheta, truethetaerror, angresthetaerror)
sigma_linearity.GetXaxis().SetTitle("#theta (deg)")
sigma_linearity.GetYaxis().SetTitle("#sigma (mrad)")
sigma_linearity.SetName("sigma_linearity")
sigma_linearity.Write()