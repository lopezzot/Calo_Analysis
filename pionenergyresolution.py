import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
import calibration
#import calibration2 as calibration

def recenergy(name):
	outputfile = "PionEnergyRes"+str(name)
	displayfile = TFile(outputfile+".root","RECREATE")
	
	if name == "FTFPBERT":
		chi = 0.3
	if name == "FTFPBERTTRV":
		chi = 0.3
	if name == "QGSPBERT":
		chi = 0.3
	if name == "QBBC":
		chi = 0.3
	
	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	Energy = array('d')
	energyfractionscin = array('d')
	energyfractioncher = array('d')
	energyfraction = array('d')
	resolutionscin = array('d')
	resolutioncher = array('d')
	resolution = array('d')

	energies = array('d',[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150])
	energies = array('d', [10,20,40,60, 80,90, 100,130, 150])
	energies = array('d', [10,30,50,70,90,100])
	sqrtenergies = array('d',[1/(x**0.5) for x in energies])
	scin_sqrtenergies = array('d')
	cher_sqrtenergies = array('d')

	t = [10,30,50,70,100,140,150]
	t = [100]
	energies = array('d',t)
	#inputfiles = ["/home/software/Calo/results/energycont_2p0m/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/software/Calo/results/pionenergyscan_QGSPBICHP/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Desktop/Calo/newresults/FTFPBERTTRV/Pion_"+str(i)+"_FTFPBERTTRV_office.root" for i in t]
	#inputfiles = ["/Users/lorenzo/Desktop/ToPC/newresults/"+str(name)+"/Pion_"+str(i)+".root" for i in t]
	inputfiles = ["/home/lorenzo/Calo/results/Pion_25_3_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t]

	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_"+str(t[counter]), str(t[counter])+"_scin", 400, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_"+str(t[counter]), str(t[counter])+"_cher", 400, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_"+str(t[counter]),str(t[counter])+"_Energy", 400, 0., 200.)
		
		#Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 3000, 0., 30000)
		EnergyHist = TH1F("Energy_"+str(t[counter]),str(t[counter])+"_Energy", 400, 0., 200.)
		LeakageHist = TH1F("Leak_"+str(t[counter]),str(t[counter])+"_Leak", 200, 0., 100.)
		NeutrinoLeakageHist = TH1F("NeutrinoNeutrinoLeak_"+str(t[counter]),str(t[counter])+"_Leak", 200, 0., 100.)
		TotalLeakageHist = TH1F("TotalLeak_"+str(t[counter]),str(t[counter])+"_Leak", 200, 0., 100.)
		ChiHist = TH1F("Chi_"+str(t[counter]),str(t[counter])+"_Chi", 40, 0., 2.)
		scatterplot = TH2F("scatterplot_"+str(t[counter]), str(t[counter]), int(400), 0., 200., int(400), 0., 200.)
		EnergyContHist = TH1F("EnergyCont_"+str(t[counter]),str(t[counter])+"_EnergyCont", 400, 0., 200.)

		#loop over events
		for Event in range(50000):	

			tree.GetEntry(Event)	
			if Event%1000==0:
				print Event 

			#Set values of the tree
			PrimaryParticleName = tree.PrimaryParticleName # MC truth: primary particle Geant4 name
			PrimaryParticleEnergy = tree.PrimaryParticleEnergy
			EnergyTot = tree.EnergyTot # Total energy deposited in calorimeter
			Energyem = tree.Energyem # Energy deposited by the em component
			EnergyScin = tree.EnergyScin # Energy deposited in Scin fibers (not Birk corrected)
			EnergyCher = tree.EnergyCher # Energy deposited in Cher fibers (not Birk corrected)
			NofCherenkovDetected = tree.NofCherenkovDetected # Total Cher p.e. detected
			BarrelR_VectorSignals = tree.VectorSignalsR  # Vector of energy deposited in Scin fibers (Birk corrected)
			BarrelL_VectorSignals = tree.VectorSignalsL  # Vector of energy deposited in Scin fibers (Birk corrected)
			BarrelR_VectorSignalsCher = tree.VectorSignalsCherR # Vector of Cher p.e. detected in Cher fibers
			BarrelL_VectorSignalsCher = tree.VectorSignalsCherL 	
			VectorR = tree.VectorR
			VectorL = tree.VectorL
			Leak = tree.leakage
			NeutrinoLeak = tree.neutrinoleakage
			
			#totalsignalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
			#Signalscinhist.Fill(totalsignalscin)

			energytot = (sum(VectorR)+sum(VectorL))/1000
			EnergyHist.Fill(energytot)
			LeakageHist.Fill(Leak/1000.)
			NeutrinoLeakageHist.Fill(NeutrinoLeak/1000.)
			TotalLeakageHist.Fill(Leak/1000.+NeutrinoLeak/1000.)


			if (Leak/1000.+NeutrinoLeak/1000.)<3.0:
				#apply calibrations
				Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
				Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
				Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
				Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)
				#end of calibrations	
				energyscin = sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals)
				energycher = sum(Calib_BarrelR_VectorSignalsCher)+sum(Calib_BarrelL_VectorSignalsCher)	
				e_c = float(t[counter])-(Leak/1000.+NeutrinoLeak/1000.)
				EnergyContHist.Fill(e_c)
				ScinEnergyHist.Fill(energyscin)
				CherEnergyHist.Fill(energycher)
				scatterplot.Fill(energyscin, energycher)
				chi = 0.29
				newchi = (energyscin-e_c)/(energycher-e_c)
				ChiHist.Fill(newchi)			
				RecEnergyHist.Fill((energyscin - chi*energycher)/(1.- chi))	
				
		print energies[counter], ScinEnergyHist.GetMean(), CherEnergyHist.GetMean(), RecEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		#ScinEnergyHist.Fit("gaus")
		#CherEnergyHist.Fit("gaus")
		RecEnergyHist.Fit("gaus")
		RecEnergyHist.Write()
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		#Signalscinhist.Write()
		EnergyHist.Write()
		EnergyContHist.Write()
		e_cont = EnergyContHist.GetMean()
		scatterplot.Write()
		LeakageHist.Write()
		NeutrinoLeakageHist.Write()
		TotalLeakageHist.Write()
		ChiHist.Write()
		#scin_sqrtenergies.append(1./(ScinEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		#cher_sqrtenergies.append(1./(CherEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		MeanEnergyScin.append(ScinEnergyHist.GetMean())
		MeanEnergyCher.append(CherEnergyHist.GetMean())
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/e_cont)
		energyfractionscin.append(ScinEnergyHist.GetMean()/e_cont)
		energyfractioncher.append(CherEnergyHist.GetMean()/e_cont)
		energyfraction.append(RecEnergyHist.GetFunction("gaus").GetParameter(1)/e_cont)
		resolutionscin.append(ScinEnergyHist.GetRMS()/e_cont)
		resolutioncher.append(CherEnergyHist.GetRMS()/e_cont)

	LinearityGraph = TGraph(len(energies), energies, energyfraction)
	LinearityGraph.SetName("LinearityGraph")
	LinearityGraph.Write()
	LinearityGraphScin = TGraph(len(energies), energies, energyfractionscin)
	LinearityGraphCher = TGraph(len(energies), energies, energyfractioncher)
	LinearityGraphCher.SetName("LinearityGraphCher")
	LinearityGraphCher.Write()
	LinearityGraphScin.SetName("LinearityGraphScin")
	LinearityGraphScin.Write()

	ResolutionGraphScin = TGraph(len(energies), energies, resolutionscin)
	func = TF1("func", "[0]/(x**0.5)+[1]", 10., 150.)
	ResolutionGraphCher = TGraph(len(energies), energies, resolutioncher)
	#ResolutionGraphScin.Fit("func", "R")
	#ResolutionGraphCher.Fit("func", "R")
	ResolutionGraphScin.SetName("ResolutionGraphScin")
	ResolutionGraphScin.Write()
	ResolutionGraphCher.SetName("ResolutionGraphCher")
	ResolutionGraphCher.Write()
	ResolutionGraph = TGraph(len(energies), energies, resolution)
	ResolutionGraph.Fit("func", "R")
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.Write()
	
	
	Resolutions = TMultiGraph()
	Resolutions.Add(ResolutionGraphScin)
	Resolutions.Add(ResolutionGraphCher)
	Resolutions.Add(ResolutionGraph)
	Resolutions.SetName("EMResolutions")
	Resolutions.Write()

	Linearities = TMultiGraph()
	Linearities.Add(LinearityGraph)
	Linearities.Add(LinearityGraphScin)
	Linearities.Add(LinearityGraphCher)
	Linearities.SetName("Linearities")
	Linearities.Write()
	
#names = ["FTFPBERT"]#, "FTFPBERT", "QGSPBERT", "QBBC"]
#names = ["FTFPBERTTRV", "QGSPBERT", "QBBC"]
name = raw_input("physics list: ")
names = []
names.append(name)

for name in names:
	recenergy(name)
