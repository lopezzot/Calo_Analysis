import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors
import glob
from array import array
import os
import newmap
import numpy as np

path = str("/home/software/Calo/results/Tower_1.root") #file to get energy in scin fibers and cherenkov photons

def getsignals():
	inputfile = TFile(path)
	tree = TTree()
	inputfile.GetObject("B4", tree)

	outputfile = raw_input("Insert root output file: ")
	file = TFile(outputfile,"RECREATE")

	SignalScin = TH1F("SignalScin", "SignalScin", 100, 0., 30000.)
	SignalCher = TH1F("SignalCher", "SignalCher", 100, 0., 50000.)
	Energydep = TH1F("EnergyDep", "Energydep", 10000, 0., 100000.)	
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
		BarrelL_VectorSignalsCher = tree.VectorSignalsCherL # Vecotr of Cher p.e. detected in Cher fibers	
		VectorR = tree.VectorR # Energy deposted
		VectorL = tree.VectorL # Energy deposited	
		
		S = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
		C = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)
		E = sum(VectorR)+sum(VectorL)
		print S
		SignalScin.Fill(S)
		SignalCher.Fill(C)
		Energydep.Fill(E)
	
	SignalScin.Write()
	SignalCher.Write()
	Energydep.Write()

	file.Close()

getsignals()