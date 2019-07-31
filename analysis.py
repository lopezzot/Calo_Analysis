import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors
import glob
from array import array
import os
import newmap
import numpy as np

machine = raw_input("On what machine running? (mac, linux, office) ")
if machine == "mac":
		path = str("/Users/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/Users/lorenzo/cernbox/work/Git-to-Mac/IDEA_Calorimeter_Union_data/")
if machine == "linux":
		path = str("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/home/lorenzo/cernbox/work/Git-to-Mac/IDEA_Calorimeter_Union_data/BarrelR/")
if machine == "office":
		datapath = str("/home/software/Calo/results/tower_calibration/")

def eventdisplay(inputfile, outputfile, histoname):
	#inputfile = raw_input("Insert root file: ")
	inputfile = TFile(datapath+inputfile)
	tree = TTree()
	inputfile.GetObject("B4", tree)

	#outputfile = raw_input("Insert root output file: ")
	displayfile = TFile(outputfile+".root","RECREATE")

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
		
		ROOTHistograms.create_eventdisplay_scin(PrimaryParticleName, BarrelR_VectorSignals, BarrelL_VectorSignals, histoname) 

def towercalibration():
	outputfile = "TowersRight"
	displayfile = TFile(outputfile+".root","RECREATE")

	MeanScin = array('d')
	MeanCher = array('d')
	RMSScin = array('d')
	RMSCher = array('d')
	Tower = array('d')
	EnergyTower = array('d')
	Zeros = array('d')
	errorsscin = array('d')
	errorscher = array('d')

	inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinHist = TH1F("scin_", str(counter+1)+"_scin", 100, 0., 10.)
		CherHist = TH1F("cher_", str(counter+1)+"_cher", 100, 0., 50.)	
		EnergyTowerHist = TH1F("Energy_",str(counter+1)+"_Energy", 200, 0., 100000.)
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

			signalscin = max(BarrelR_VectorSignals)
			signalcher = max(BarrelR_VectorSignalsCher)
			energytower = max(VectorR)
			
			ScinHist.Fill(energytower/signalscin)	
			CherHist.Fill(energytower/signalcher)	
			EnergyTowerHist.Fill(energytower)

			if Event < 1:
				print "Max found at: "+str(list(BarrelR_VectorSignals).index(signalscin))+str(list(BarrelR_VectorSignalsCher).index(signalcher))+str(list(VectorR).index(energytower))+" for file "+str(counter+1) #to check tower mostly hitten is the correct one
				displayfile.cd()
				ROOTHistograms.create_eventdisplay_scin(PrimaryParticleName, BarrelR_VectorSignals, BarrelL_VectorSignals, str(counter)) 
		
		displayfile.cd()
		ScinHist.Write()
		CherHist.Write()
		EnergyTowerHist.Write()
		ScinHist.Fit("gaus")
		MeanScin.append(ScinHist.GetMean())
		MeanCher.append(CherHist.GetMean())
		RMSScin.append(ScinHist.GetRMS())
		RMSCher.append(CherHist.GetRMS())
		EnergyTower.append(EnergyTowerHist.GetMean())
		#print ScinHist.GetMean(), CherHist.GetMean(), EnergyTowerHist.GetMean()
		#print ScinHist.GetRMS(), CherHist.GetRMS()
		Tower.append(counter+1)
		Zeros.append(0.)
		errorsscin.append(ScinHist.GetRMS()/(3000**0.5))
		errorscher.append(CherHist.GetRMS()/(3000**0.5))

	#print np.mean(MeanScin), max(MeanScin), min(MeanScin)
	#print np.mean(MeanCher), max(MeanCher), min(MeanCher)
	
	n = len(MeanCher)
	RMSGraphScin = TGraph(n, Tower, RMSScin)
	RMSGraphCher = TGraph(n, Tower, RMSCher)
	MeanGraphScin = TGraphErrors(n, Tower, MeanScin, Zeros, errorsscin)
	MeanGraphCher = TGraphErrors(n, Tower, MeanCher, Zeros, errorscher)
	x = array('d', (0., 90., 90., 0.))
	y = array('d', (np.mean(MeanScin)-0.015*np.mean(MeanScin), np.mean(MeanScin)-0.015*np.mean(MeanScin), np.mean(MeanScin)+0.015*np.mean(MeanScin), np.mean(MeanScin)+0.015*np.mean(MeanScin)))
	Fillgraph = TGraph(4, x, y )
	Fillgraph.Write()
	x2 = array('d', (0., 90., 90., 0.))
	y2 = array('d', (np.mean(MeanCher)-0.01*np.mean(MeanCher), np.mean(MeanCher)-0.01*np.mean(MeanCher), np.mean(MeanCher)+0.01*np.mean(MeanCher), np.mean(MeanCher)+0.01*np.mean(MeanCher)))
	Fillgraph2 = TGraph(4, x2, y2 )
	Fillgraph2.Write()
	MeanGraphCher.Write()
	MeanGraphScin.Write()
	RMSGraphCher.Write()
	RMSGraphScin.Write()

towercalibration()