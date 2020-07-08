import os
import numpy as np
#import matplotlib.pyplot as plt
import time
from ROOT import TH2F, TFile, TH1F, TGraph2D, TLine, gStyle, TGraph, TVector3
import circle
import createangularesolutionplot
import math
'''
def cart2sph(x,y,z):
	phi = np.arctan2(y,x)
	theta = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return theta, phi, r
'''
def cart2sph(x,y,z):
	if y >= 0.:
		phi = np.arccos(x/np.sqrt(x**2+y**2))
	else:
		phi = 2*math.pi-np.arccos(x/np.sqrt(x**2+y**2))
	#phi = np.arctan2(y,x)
	theta = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return theta, phi, r

particle = raw_input("particle: ")

energies = [10,20,30,40,50,60,70,80,90,100,130,140,150]
energies = [10]
outputfile = TFile(particle+"Angle.root", "RECREATE")
for e in energies:

	file = "/home/software/Calo/results/energy_angular_pion/Pion_"+str(e)+".txt"	
	file = "/home/software/Calo/results/Electron_ang_res_1_1/Electron_"+str(e)+".txt"
	file = "/home/software/Calo/NewResults/AngleRes_"+particle+"/Energy_"+str(e)+"_ThetaPhi_1.0_1.0.txt"
	#file = "/home/software/Calo/NewResults/AngleRes_"+particle+"_2/Energy_"+str(e)+".txt"
	file = "/home/software/Calo/NewResults/AngleRes_images/Energy_40_"+str(particle)+"_ThetaPhi_1.0_180.0.txt"
	file = "/home/software/Calo/NewResults/AngleRes_images/Energy_jet_90.txt"
	if e == energies[0]:
		energies = []
		angrestheta = []
		angresthetaerror = []
		angresphi = []
		angresphierror = []
		angrestheta_c = []
		angrestheta_c_error = []
		angrestheta_s = []
		angrestheta_s_error = []
		angresphi_c = []
		angresphi_c_error = []
		angresphi_s = []
		angresphi_s_error = []
		fiber_s = []
		fiber_s_error = []
		fiber_c = []
		fiber_c_error = []
		fiber_s_1suppression = []
		fiber_c_1suppression = []
		fiber_s_2suppression = []
		fiber_c_2suppression = []

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

	percentages_array_s = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
	percentages_array_c = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])


	FiberSHist = TH1F("FiberSHist_"+str(e), "FiberS_"+str(e),10000,0.,10000.)
	FiberCHist = TH1F("FiberCHist_"+str(e), "FiberC_"+str(e),10000,0.,10000.)
	FiberSHist_1suppression = TH1F("FiberSHist_"+str(e)+"_1sup", "FiberS_"+str(e),10000,0.,10000.)
	FiberCHist_1suppression = TH1F("FiberCHist_"+str(e)+"_1sup", "FiberC_"+str(e),10000,0.,10000.)
	FiberSHist_2suppression = TH1F("FiberSHist_"+str(e)+"_2sup", "FiberS_"+str(e),10000,0.,10000.)
	FiberCHist_2suppression = TH1F("FiberCHist_"+str(e)+"_2sup", "FiberC_"+str(e),10000,0.,10000.)
	for i in range (0,int(max(EvtID[1:]))+1):
		# S have Flag == 1 C == 0
		S = np.where((Flag==1.) & (EvtID==i))
		X_S=X[S]
		Y_S=Y[S]
		Z_S=Z[S]
		theta_S=[]
		phi_S=[]
		E_s = Edep[S]
		FiberSHist.Fill(len(X_S))
		FiberSHist_1suppression.Fill(len([x for x in E_s if x>1.0]))
		FiberSHist_2suppression.Fill(len([x for x in E_s if x>2.0]))
		C = np.where((Flag==0.) & (EvtID==i))
		X_C=X[C]
		Y_C=Y[C]
		Z_C=Z[C]
		theta_C=[]
		phi_C=[]
		E_c=Edep[C]	
		FiberCHist.Fill(len(X_C))
		FiberCHist_1suppression.Fill(len([x for x in E_c if x>1.0]))
		FiberCHist_2suppression.Fill(len([x for x in E_c if x>2.0]))

		for counter, j in enumerate(X_S):
			(Stheta, Sphi, Sr)=cart2sph(float(X_S[counter]),float(Y_S[counter]),float(Z_S[counter]))
			theta_S.append(Stheta)
			phi_S.append(Sphi)	
				
		for counter, j in enumerate(X_C):
			(Ctheta, Cphi, Cr)=cart2sph(float(X_C[counter]),float(Y_C[counter]),float(Z_C[counter]))
			theta_C.append(Ctheta)
			phi_C.append(Cphi)		

		#ScinPlot = TH2F("Scin_"+str(e), "Scin_"+str(e), int(1000), float(np.mean(theta_S)-0.5), float(np.mean(theta_S)+0.5), int(1000), float(np.mean(phi_S)-0.5),float(np.mean(phi_S)+0.5))
		ScinPlot = TH2F("Scin_"+str(e)+"_"+str(i), "Scin_"+str(e), int(1000), float(np.min(theta_S)), float(np.max(theta_S)), int(1000), float(np.min(phi_S)),float(np.max(phi_S)))
		CherPlot = TH2F("Cher_"+str(e)+"_"+str(i), "Cher_"+str(e), int(1000), float(np.min(theta_S)), float(np.max(theta_S)), int(1000), float(np.min(phi_S)),float(np.max(phi_S)))
		
		n = len(theta_S)
		array_theta_S = np.array(theta_S, 'd')
		array_phi_S = np.array(phi_S, 'd')
		array_E_s = np.array(E_s, 'd')
		array_E_s_ordered, array_theta_S_ordered, array_phi_S_ordered = zip(*sorted(zip(array_E_s, array_theta_S, array_phi_S)))
		array_theta_S_ordered = [x*180./math.pi+90. for x in array_theta_S_ordered]
		array_phi_S_ordered = [x*180./math.pi for x in array_phi_S_ordered]
		array_theta_S_ordered = array_theta_S_ordered + [0.,0.,180.,180.]
		array_phi_S_ordered = array_phi_S_ordered + [0.,360.,0.,360.]
		array_E_s_ordered = [x for x in array_E_s_ordered] + [-100.,-100.,-100.,-100.]
		ScinGraph = TGraph2D(n+4, np.array(array_theta_S_ordered), np.array(array_phi_S_ordered), np.array(array_E_s_ordered))
		ScinGraph.SetMinimum(0.0)
		ScinGraph.SetName("ScinGraph_"+str(e)+"_"+str(i))

		n = len(theta_C)
		array_theta_C = np.array(theta_C, 'd')
		array_phi_C = np.array(phi_C, 'd')
		array_E_c = np.array(E_c, 'd')
		array_E_c_ordered, array_theta_C_ordered, array_phi_C_ordered = zip(*sorted(zip(array_E_c, array_theta_C, array_phi_C)))
		array_theta_C_ordered = [x*180./math.pi+90. for x in array_theta_C_ordered]
		array_phi_C_ordered = [x*180./math.pi for x in array_phi_C_ordered]
		array_theta_C_ordered = array_theta_C_ordered + [0.,0.,180.,180.]
		array_phi_C_ordered = array_phi_C_ordered + [0.,360.,0.,360.]
		array_E_c_ordered = [x for x in array_E_c_ordered] + [-100.,-100.,-100.,-100.]
		CherGraph = TGraph2D(n+4, np.array(array_theta_C_ordered), np.array(array_phi_C_ordered), np.array(array_E_c_ordered))
		CherGraph.SetMinimum(0.0)
		CherGraph.SetName("CherGraph_"+str(e)+"_"+str(i))

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
			#if -0.04<theta_S[counter]<0.08:
			MeanTheta += entry*theta_S[counter]
			sumtheta += entry
			#if -0.04<phi_S[counter]<0.08:
			MeanPhi += entry*phi_S[counter]
			sumphi += entry
			ScinPlot.Fill(theta_S[counter], phi_S[counter], entry)	

		for counter, entry in enumerate(E_c):
			#if -0.04<theta_C[counter]<0.08:
			MeanThetaC += entry*theta_C[counter]
			sumthetaC += entry
			#if -0.04<phi_C[counter]<0.08:
			MeanPhiC += entry*phi_C[counter]
			sumphiC += entry
			CherPlot.Fill(theta_C[counter], phi_C[counter], entry)	

		for counter, entry in enumerate(E_s):
			#if -0.04<theta_S[counter]<0.08:
			MeanThetaCS += entry*theta_S[counter]
			sumthetaCS += entry
			#if -0.04<phi_S[counter]<0.08:
			MeanPhiCS += entry*phi_S[counter]
			sumphiCS += entry
		for counter, entry in enumerate(E_c):
			#if -0.04<theta_C[counter]<0.08:
			MeanThetaCS += entry*theta_C[counter]
			sumthetaCS += entry
			#if -0.04<phi_C[counter]<0.08:
			MeanPhiCS += entry*phi_C[counter]
			sumphiCS += entry
			
		MeanTheta = MeanTheta/sumtheta
		MeanPhi = MeanPhi/sumphi

		ScinPlot.GetXaxis().SetRangeUser(-0.04,0.08)
		ScinPlot.GetYaxis().SetRangeUser(-0.04,0.08)
		CherPlot.GetXaxis().SetRangeUser(-0.04,0.08)
		CherPlot.GetYaxis().SetRangeUser(-0.04,0.08)		
		MeanTheta2 = ScinPlot.GetMean(1)
		#print MeanTheta, MeanTheta2, MeanTheta-MeanTheta2
		ThetaHist.Fill(MeanTheta)
		PhiHist.Fill(MeanPhi)
		percentages_s = circle.computerradi(E_s, theta_S, phi_S)
		
		MeanThetaC = MeanThetaC/sumthetaC
		MeanPhiC = MeanPhiC/sumphiC
		ThetaHistC.Fill(MeanThetaC)
		PhiHistC.Fill(MeanPhiC)	
		percentages_c = circle.computerradi(E_c, theta_C, phi_C)

		MeanThetaCS = MeanThetaCS/sumthetaCS
		MeanPhiCS = MeanPhiCS/sumphiCS
		ThetaHistCS.Fill(MeanThetaCS)
		PhiHistCS.Fill(MeanPhiCS)	

		percentages_array_s = percentages_array_s + np.array(percentages_s)
		percentages_array_c = percentages_array_c + np.array(percentages_c) 
		if i<100:
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

	FiberSHist.Write()
	FiberCHist.Write()
	fiber_s.append(FiberSHist.GetMean())
	fiber_s_error.append(FiberSHist.GetRMS()/(3000**0.5))
	fiber_c.append(FiberCHist.GetMean())
	fiber_c_error.append(FiberCHist.GetRMS()/(3000**0.5))
	fiber_s_1suppression.append(FiberSHist_1suppression.GetMean())
	fiber_s_2suppression.append(FiberSHist_2suppression.GetMean())
	fiber_c_1suppression.append(FiberCHist_1suppression.GetMean())
	fiber_c_2suppression.append(FiberCHist_2suppression.GetMean())

	percentages_array_s = percentages_array_s/3000
	percentages_array_c = percentages_array_c/3000
	print percentages_array_s
	radii = np.array([0.001, 0.003, 0.005, 0.01, 0.05, 0.1, 0.15], 'd')
	p_s = TGraph(len(radii), radii, percentages_array_s)
	p_s.SetName("Radius_s_"+str(e))
	p_s.SetTitle("radius_s_"+str(e))
	p_s.Write()

	p_c = TGraph(len(radii), radii, percentages_array_c)
	p_c.SetName("Radius_c_"+str(e))
	p_c.SetTitle("radius_c_"+str(e))
	p_c.Write()

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
	angresthetaerror.append(ThetaHistCS.GetFunction("gaus").GetParError(2)*1000)
	angresphi.append(PhiHistCS.GetFunction("gaus").GetParameter(2)*1000)
	angresphierror.append(PhiHistCS.GetFunction("gaus").GetParError(2)*1000)	
	angrestheta_s.append(ThetaHist.GetFunction("gaus").GetParameter(2)*1000)
	angrestheta_s_error.append(ThetaHist.GetFunction("gaus").GetParError(2)*1000)
	angresphi_s.append(PhiHist.GetFunction("gaus").GetParameter(2)*1000)
	angresphi_s_error.append(PhiHist.GetFunction("gaus").GetParError(2)*1000)
	angrestheta_c.append(ThetaHistC.GetFunction("gaus").GetParameter(2)*1000)
	angrestheta_c_error.append(ThetaHistC.GetFunction("gaus").GetParError(2)*1000)	
	angresphi_c.append(PhiHistC.GetFunction("gaus").GetParameter(2)*1000)
	angresphi_c_error.append(PhiHistC.GetFunction("gaus").GetParError(2)*1000)
		
	print angrestheta, angresphi, energies

gStyle.SetOptStat(111)
ThetaGraph, PhiGraph = createangularesolutionplot.angresplot(energies, angrestheta, angresphi, angresthetaerror, angresphierror)
ThetaGraph_s, PhiGraph_s = createangularesolutionplot.angresplot(energies, angrestheta_s, angresphi_s, angrestheta_s_error, angresphi_s_error)
ThetaGraph_c, PhiGraph_c = createangularesolutionplot.angresplot(energies, angrestheta_c, angresphi_c, angrestheta_c_error, angresphi_c_error)
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

Fiber_S_graph = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_s,'d'))
Fiber_S_graph.SetName("Fiber_S")
Fiber_S_graph.Write()
Fiber_C_graph = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_c,'d'))
Fiber_C_graph.SetName("Fiber_C")
Fiber_C_graph.Write()
Fiber_S_graph_1suppression = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_s_1suppression,'d'))
Fiber_S_graph_1suppression.SetName("fiber_s_1suppression")
Fiber_S_graph_1suppression.Write()
Fiber_C_graph_1suppression = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_c_1suppression,'d'))
Fiber_C_graph_1suppression.SetName("fiber_c_1suppression")
Fiber_C_graph_1suppression.Write()
Fiber_S_graph_2suppression = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_s_2suppression,'d'))
Fiber_S_graph_2suppression.SetName("fiber_s_2suppression")
Fiber_S_graph_2suppression.Write()
Fiber_C_graph_2suppression = TGraph(len(energies), np.array(energies,'d'), np.array(fiber_c_2suppression,'d'))
Fiber_C_graph_2suppression.SetName("fiber_c_2suppression")
Fiber_C_graph_2suppression.Write()



