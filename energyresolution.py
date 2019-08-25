import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1
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
		datapath = str("/home/software/Calo/results/NewTowerScan4/")

def recenergy():
	outputfile = "EMEnergyRes"
	displayfile = TFile(outputfile+".root","RECREATE")

	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	Energy = array('d')
	
	inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	#inputfiles = ["/home/software/Calo/results/NewTowerScan4/Barrel_"+str(i)+".root" for i in range(1,76)]
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 200, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 200, 0., 200.)	
		EnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 200, 0., 200.)
		
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

		print ScinEnergyHist.GetMean(), CherEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		ScinEnergyHist.Fit("gaus")
		CherEnergyHist.Fit("gaus")
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		Signalscinhist.Write()
		EnergyHist.Write()
		MeanEnergyScin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		MeanEnergyCher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1))

recenergy()