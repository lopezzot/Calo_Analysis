import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
import calibration

machine = raw_input("On what machine running? (mac, linux, office) ")
if machine == "mac":
		path = str("/Users/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/Users/lorenzo/cernbox/work/Git-to-Mac/IDEA_Calorimeter_Union_data/")
if machine == "linux":
		path = str("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/home/lorenzo/Desktop/Calo/results/newenergyscan3/")
if machine == "office":
		datapath = str("/home/software/Calo/results/pionenergyscan/")

def recenergy():
	outputfile = "PionEnergyRes"
	displayfile = TFile(outputfile+".root","RECREATE")

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
	energies = array('d', [100])
	sqrtenergies = array('d',[1/(x**0.5) for x in energies])
	scin_sqrtenergies = array('d')
	cher_sqrtenergies = array('d')
	#inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	t = [10,20, 40, 60, 80,90, 100, 130, 150]
	t = [100]
	#inputfiles = ["/home/software/Calo/results/energycont_2p0m/Pion_"+str(i)+".root" for i in t]

	inputfiles = ["/home/software/Calo/results/pionenergyscan_QGSPBICHP/Pion_"+str(i)+".root" for i in t]
	inputfiles = ["/home/lorenzo/Desktop/Calo/results/pionenergyscan_FTFPBERTTRV/Pion_"+str(i)+".root" for i in t]

	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 500, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 500, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 500, 0., 200.)
		
		#Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 3000, 0., 30000)
		EnergyHist = TH1F("Energy_",str(counter+1)+"_Energy", 500, 0., 200.)
		
		scatterplot = TH2F("scatterplot_", str(counter+1), int(800), 0., 200., int(800), 0., 200.)

		#loop over events
		for Event in range(int(tree.GetEntries())):	

			tree.GetEntry(Event)	

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
			
			#totalsignalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
			#Signalscinhist.Fill(totalsignalscin)

			energytot = (sum(VectorR)+sum(VectorL))/1000
			EnergyHist.Fill(energytot)
			
			#apply calibrations
			Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
			Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
			Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
			Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)
			#end of calibrations

			energyscin = sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals)
			energycher = sum(Calib_BarrelR_VectorSignalsCher)+sum(Calib_BarrelL_VectorSignalsCher)

			ScinEnergyHist.Fill(energyscin)
			#sigmascin = 0.15*(energyscin**0.5)+0.012*energyscin
			CherEnergyHist.Fill(energycher)
			#sigmacher = 0.18*(energycher**0.5)+0.0045*energycher
			chi = 0.29
			RecEnergyHist.Fill((energyscin - chi*energycher)/(1-chi))

			scatterplot.Fill(energyscin, energycher)
		
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
		scatterplot.Write()
		#scin_sqrtenergies.append(1./(ScinEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		#cher_sqrtenergies.append(1./(CherEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		MeanEnergyScin.append(ScinEnergyHist.GetMean())
		MeanEnergyCher.append(CherEnergyHist.GetMean())
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))
		energyfractionscin.append(ScinEnergyHist.GetMean()/energies[counter])
		energyfractioncher.append(CherEnergyHist.GetMean()/energies[counter])
		energyfraction.append(RecEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		#resolutionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		#resolutioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(2)/CherEnergyHist.GetFunction("gaus").GetParameter(1))

	LinearityGraph = TGraph(len(energies), energies, energyfraction)
	LinearityGraph.SetName("LinearityGraph")
	LinearityGraph.Write()
	LinearityGraphScin = TGraph(len(energies), energies, energyfractionscin)
	LinearityGraphCher = TGraph(len(energies), energies, energyfractioncher)
	LinearityGraphCher.SetName("LinearityGraphCher")
	LinearityGraphCher.Write()
	LinearityGraphScin.SetName("LinearityGraphScin")
	LinearityGraphScin.Write()

	#ResolutionGraphScin = TGraph(len(energies), scin_sqrtenergies, resolutionscin)
	func = TF1("func", "[0]/(x**0.5)+[1]", 10., 150.)
	#ResolutionGraphCher = TGraph(len(energies), cher_sqrtenergies, resolutioncher)
	#ResolutionGraphScin.Fit("func", "R")
	#ResolutionGraphCher.Fit("func", "R")
	#ResolutionGraphScin.SetName("ResolutionGraphScin")
	#ResolutionGraphScin.Write()
	#ResolutionGraphCher.SetName("ResolutionGraphCher")
	#ResolutionGraphCher.Write()
	ResolutionGraph = TGraph(len(energies), energies, resolution)
	ResolutionGraph.Fit("func", "R")
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.Write()
	'''
	rd52copper = array('d', [0.04478505426185217, 0.027392527130926082, 0.02420093893609386, 0.02229837387624884, 0.020999999999999998])

	rd52graph = TGraph(len(energies), sqrtenergies, rd52copper)
	rd52graph.SetName("rd52resolution")
	rd52graph.Write()


	EMResolutions = TMultiGraph()
	EMResolutions.Add(ResolutionGraphScin)
	EMResolutions.Add(ResolutionGraphCher)
	EMResolutions.Add(ResolutionGraph)
	EMResolutions.Add(rd52graph)
	EMResolutions.SetName("EMResolutions")
	EMResolutions.Write()

	Linearities = TMultiGraph()
	Linearities.Add(LinearityGraph)
	Linearities.Add(LinearityGraphScin)
	Linearities.Add(LinearityGraphCher)
	Linearities.SetName("Linearities")
	Linearities.Write()
	'''

recenergy()
