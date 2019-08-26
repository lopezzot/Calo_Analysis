import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1, TMultiGraph
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
		datapath = str("/home/software/Calo/results/newenergyscan3/")

def recenergy():
	outputfile = "EMEnergyRes"
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
	sqrtenergies = array('d',[1/(x**0.5) for x in energies])

	inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	#inputfiles = ["/home/software/Calo/results/NewTowerScan4/Barrel_"+str(i)+".root" for i in range(1,76)]
	for counter, inputfile in enumerate(inputfiles[1:]):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 200, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 200, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 200, 0., 200.)
		
		Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 3000, 0., 30000)
		EnergyHist = TH1F("Energy_",str(counter+1)+"_Energy", 200, 0., 200.)
		
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
			
			totalsignalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
			Signalscinhist.Fill(totalsignalscin)

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
			CherEnergyHist.Fill(energycher)
			RecEnergyHist.Fill((energyscin+energycher)/2)

		print ScinEnergyHist.GetMean(), CherEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		ScinEnergyHist.Fit("gaus")
		CherEnergyHist.Fit("gaus")
		RecEnergyHist.Fit("gaus")
		RecEnergyHist.Write()
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		Signalscinhist.Write()
		EnergyHist.Write()
		MeanEnergyScin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		MeanEnergyCher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1))
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))
		energyfractionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		energyfractioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		energyfraction.append(RecEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		resolutionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		resolutioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(2)/CherEnergyHist.GetFunction("gaus").GetParameter(1))

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

	func = TF1("func", "[0]/(x**0.5)+[1]", 10, 100)
	ResolutionGraphCher = TGraph(len(energies), energies, resolutioncher)
	ResolutionGraphScin.Fit("func", "R")
	ResolutionGraphCher.Fit("func", "R")
	ResolutionGraphScin.SetName("ResolutionGraphScin")
	ResolutionGraphScin.Write()
	ResolutionGraphCher.SetName("ResolutionGraphCher")
	ResolutionGraphCher.Write()
	ResolutionGraph = TGraph(len(energies), energies, resolution)
	ResolutionGraph.Fit("func", "R")
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.Write()

	EMResolutions = TMultiGraph()
	EMResolutions.Add(ResolutionGraphScin)
	EMResolutions.Add(ResolutionGraphCher)
	EMResolutions.Add(ResolutionGraph)
	EMResolutions.SetName("EMResolutions")
	EMResolutions.Write()

	Linearities = TMultiGraph()
	Linearities.Add(LinearityGraph)
	Linearities.Add(LinearityGraphScin)
	Linearities.Add(LinearityGraphCher)
	Linearities.SetName("Linearities")
	Linearities.Write()

def energylinearity():
	outputfile = "EMLinearityEnergyRes"
	displayfile = TFile(outputfile+".root","RECREATE")

	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	MeanEnergy = array('d')
	resolutionscin = array('d')
	resolutioncher = array('d')
	resolution = array('d')

	inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	#inputfiles = ["/home/software/Calo/results/NewTowerScan4/Barrel_"+str(i)+".root" for i in range(1,76)]
	for counter, inputfile in enumerate(inputfiles[1:]):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 200, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 200, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 200, 0., 200.)
		
		energy = 40.0
		sqrtenergy = 1/(40.0**0.5)
		towers = array('d', range(1,76))
	
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
			
			#apply calibrations
			Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
			Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
			Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
			Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)
			#end of calibrations

			energyscin = sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals)
			energycher = sum(Calib_BarrelR_VectorSignalsCher)+sum(Calib_BarrelL_VectorSignalsCher)

			ScinEnergyHist.Fill(energyscin)
			CherEnergyHist.Fill(energycher)
			RecEnergyHist.Fill((energyscin+energycher)/2)

		print ScinEnergyHist.GetMean(), CherEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		ScinEnergyHist.Fit("gaus")
		CherEnergyHist.Fit("gaus")
		RecEnergyHist.Fit("gaus")
		RecEnergyHist.Write()
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		MeanEnergyScin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1)/energy)
		MeanEnergyCher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1)/energy)
		MeanEnergy.append(EnergyHist.GetFunction("gaus").GetParameter(1)/energy)
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))
		resolutionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		resolutioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(2)/CherEnergyHist.GetFunction("gaus").GetParameter(1))

	LinearityGraph = TGraph(len(towers), towers, MeanEnergy)
	LinearityGraph.SetName("LinearityGraph")
	LinearityGraph.Write()
	LinearityGraphScin = TGraph(len(towers), towers, MeanEnergyScin)
	LinearityGraphCher = TGraph(len(towers), towers, MeanEnergyCher)
	LinearityGraphCher.SetName("LinearityGraphCher")
	LinearityGraphCher.Write()
	LinearityGraphScin.SetName("LinearityGraphScin")
	LinearityGraphScin.Write()
	ResolutionGraphScin = TGraph(len(towers), towers, resolutionscin)
	ResolutionGraphCher = TGraph(len(towers), towers, resolutioncher)
	ResolutionGraphScin.SetName("ResolutionGraphScin")
	ResolutionGraphScin.Write()
	ResolutionGraphCher.SetName("ResolutionGraphCher")
	ResolutionGraphCher.Write()
	ResolutionGraph = TGraph(len(towers), towers, resolution)
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.Write()

	EMResolutions = TMultiGraph()
	EMResolutions.Add(ResolutionGraphScin)
	EMResolutions.Add(ResolutionGraphCher)
	EMResolutions.Add(ResolutionGraph)
	EMResolutions.SetName("EMResolutions")
	EMResolutions.Write()

	Linearities = TMultiGraph()
	Linearities.Add(LinearityGraph)
	Linearities.Add(LinearityGraphScin)
	Linearities.Add(LinearityGraphCher)
	Linearities.SetName("Linearities")
	Linearities.Write()


#energylinearity()
recenergy()