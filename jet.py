import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1
import glob
from array import array
import os
import newmap
import numpy as np

def jetdisplay():
	outputfile = "Jetdisplay"
	displayfile = TFile(outputfile+".root","RECREATE")

	inputfile = "/home/software/Calo/results/hepmc_zll_FTFPBERTTRV/zll.root"
	#inputfiles = ["/home/lorenzo/Desktop/Calo/results/NewTowerScan4/Barrel_"+str(i)+".root" for i in range(1,76)]
	#inputfile = "/home/software/Calo/results/NewTowerScan4/Barrel_10.root"
	#inputfile = "/home/software/Calo/results/SliceScan/Slice_10.root"

	inputfile = TFile(inputfile)
	print "Analyzing: "+str(inputfile)+" \n"
	tree = TTree()
	inputfile.GetObject("B4", tree)	

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
		
		if Event < 10:
				displayfile.cd()
				ROOTHistograms.create_eventdisplay_scin("Jet", BarrelR_VectorSignals, BarrelL_VectorSignals, str(Event)) 
				ROOTHistograms.create_eventdisplay_cher("Jet", BarrelR_VectorSignalsCher, BarrelL_VectorSignalsCher, str(Event))


jetdisplay()