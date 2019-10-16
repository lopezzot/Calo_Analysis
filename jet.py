import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1, TLorentzVector
import glob
from array import array
import os
import newmap
import numpy as np
import fastjet
import calibration
import math

def jetdisplay():
	outputfile = "Jetdisplay"
	displayfile = TFile(outputfile+".root","RECREATE")

	inputfile = "zjj.root"
	#inputfiles = ["/home/lorenzo/Desktop/Calo/results/NewTowerScan4/Barrel_"+str(i)+".root" for i in range(1,76)]
	#inputfile = "/home/software/Calo/results/NewTowerScan4/Barrel_10.root"
	#inputfile = "/home/software/Calo/results/SliceScan/Slice_10.root"

	inputfile = TFile(inputfile)
	print "Analyzing: "+str(inputfile)+" \n"
	tree = TTree()
	inputfile.GetObject("B4", tree)	

	#loop over events
	for Event in range(int(100)):	

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
		
		Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
		Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
		Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
		Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)

		energy = float(sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals))
		
		threshold = 1.0 #(GeV)	

		if energy>50.:
			
			#event displays with signals (p.e.)
			#if Event < 1:
				#displayfile.cd()
				#ROOTHistograms.create_eventdisplay_scin("Jet", BarrelR_VectorSignals, BarrelL_VectorSignals, "signal"+str(Event)) 
				#ROOTHistograms.create_eventdisplay_cher("Jet", BarrelR_VectorSignalsCher, BarrelL_VectorSignalsCher, "signal"+str(Event))
		
			#event displays with energy (GeV)
			
			if True:
					displayfile.cd()
					ROOTHistograms.create_eventdisplay_scin("Jet_energy", Calib_BarrelR_VectorSignals, Calib_BarrelL_VectorSignals, "energy"+str(Event), threshold) 
					ROOTHistograms.create_eventdisplay_cher("Jet_energy", Calib_BarrelR_VectorSignalsCher, Calib_BarrelL_VectorSignalsCher, "energy"+str(Event), threshold)	
			
			inputparticles_scin = []

			#right part
			for towerindex in range(75*36):	
				theta, phi, eta = newmap.maptower(towerindex, "right")
				energy_scin = Calib_BarrelR_VectorSignals[towerindex]
				pt_scin = energy_scin*np.sin(theta*math.pi/180.)

				energy_cher = Calib_BarrelR_VectorSignalsCher[towerindex]
				pt_cher = energy_cher*np.sin(theta*math.pi/180.)
				
				towerscin = TLorentzVector()
				towerscin.SetPtEtaPhiM(pt_scin, eta, phi*math.pi/180., 0.)
				towercher = TLorentzVector()
				towercher.SetPtEtaPhiM(pt_cher, eta, phi*math.pi/180., 0.)	

				if energy_scin > threshold:
					#print "truth towers: "+str(energy_scin)+" "+str(eta)+" "+str(phi)
					inputparticles_scin.append(fastjet.PseudoJet(towerscin.Px(), towerscin.Py(), towerscin.Pz(), towerscin.E()))
				
			#left part
			for towerindex in range(75*36):	
				theta, phi, eta = newmap.maptower(towerindex, "left")
				energy_scin = Calib_BarrelL_VectorSignals[towerindex]
				pt_scin = energy_scin*np.sin(theta*math.pi/180.)
				
				energy_cher = Calib_BarrelL_VectorSignalsCher[towerindex]
				pt_cher = energy_cher*np.sin(theta*math.pi/180.)
				
				towerscin = TLorentzVector()
				towerscin.SetPtEtaPhiM(pt_scin, eta, phi*math.pi/180., 0.)
				towercher = TLorentzVector()
				towercher.SetPtEtaPhiM(pt_cher, eta, phi*math.pi/180., 0.)	

				if energy_scin > threshold:
					#print "truth towers: "+str(energy_scin)+" "+str(eta)+" "+str(phi)
					inputparticles_scin.append(fastjet.PseudoJet(towerscin.Px(), towerscin.Py(), towerscin.Pz(), towerscin.E()))
			
			jet_def = fastjet.JetDefinition(fastjet.kt_algorithm, 1.0)

			clust_seq = fastjet.ClusterSequence(inputparticles_scin, jet_def)	

			print "Event: "+str(Event)+" energy (GeV): "+str(energy)+" n-jets: "+str(len(clust_seq.inclusive_jets()))+" truth: "+str(len(inputparticles_scin))
			'''
			for jet in clust_seq.inclusive_jets():
				print "***********"
				print jet.e(), jet.eta(), jet.phi()
				for cjet in jet.constituents():
					print cjet.e(), cjet.eta(), cjet.phi()
			'''
			'''
			if len(clust_seq.inclusive_jets()) == 3:
				displayfile.cd()
				ROOTHistograms.create_eventdisplay_scin("Jet_energy", Calib_BarrelR_VectorSignals, Calib_BarrelL_VectorSignals, "energy"+str(Event), threshold) 
				ROOTHistograms.create_eventdisplay_cher("Jet_energy", Calib_BarrelR_VectorSignalsCher, Calib_BarrelL_VectorSignalsCher, "energy"+str(Event), threshold)	
				print "*****"
				for jet in clust_seq.inclusive_jets():
					print jet.e(), jet.eta(), jet.phi()
			'''
jetdisplay()