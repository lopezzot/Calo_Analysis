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

def mergejet(jet1_scin, jet2_scin, jet1_cher, jet2_cher):
	
	print "merging clusters"

	deltar1 = jet1_scin.delta_R(jet1_cher)
	deltar2 = jet1_scin.delta_R(jet2_cher)

	c = 0.34 #chi factor

	if deltar1 < deltar2:
		jet1Px = (jet1_scin.px()-c*jet1_cher.px())/(1-c)
		jet1Py = (jet1_scin.py()-c*jet1_cher.py())/(1-c)
		jet1Pz = (jet1_scin.pz()-c*jet1_cher.pz())/(1-c)
		jet1E = (jet1_scin.e()-c*jet1_cher.e())/(1.-c)
		jet1 = fastjet.PseudoJet(jet1Px, jet1Py, jet1Pz, jet1E)

		jet2Px = (jet2_scin.px()-c*jet2_cher.px())/(1-c)
		jet2Py = (jet2_scin.py()-c*jet2_cher.py())/(1-c)
		jet2Pz = (jet2_scin.pz()-c*jet2_cher.pz())/(1-c)
		jet2E = (jet2_scin.e()-c*jet2_cher.e())/(1.-c)
		jet2 = fastjet.PseudoJet(jet2Px, jet2Py, jet2Pz, jet2E)

		
	else:
		jet1Px = (jet1_scin.px()-c*jet2_cher.px())/(1-c)
		jet1Py = (jet1_scin.py()-c*jet2_cher.py())/(1-c)
		jet1Pz = (jet1_scin.pz()-c*jet2_cher.pz())/(1-c)
		jet1E = (jet1_scin.e()-c*jet2_cher.e())/(1.-c)
		jet1 = fastjet.PseudoJet(jet1Px, jet1Py, jet1Pz, jet1E)

		jet2Px = (jet2_scin.px()-c*jet1_cher.px())/(1-c)
		jet2Py = (jet2_scin.py()-c*jet1_cher.py())/(1-c)
		jet2Pz = (jet2_scin.pz()-c*jet1_cher.pz())/(1-c)
		jet2E = (jet2_scin.e()-c*jet1_cher.e())/(1.-c)
		jet2 = fastjet.PseudoJet(jet2Px, jet2Py, jet2Pz, jet2E)

	return jet1, jet2
#used to calculate w mass including muon (i.e. muon is not subtracted from truth)
def jetdisplay():
	outputfile = "Jetdisplay"
	displayfile = TFile(outputfile+".root","RECREATE")

	inputfile = "wwlj1k.root"

	inputfile = TFile(inputfile)
	print "Analyzing: "+str(inputfile)+" \n"
	tree = TTree()
	inputfile.GetObject("B4", tree)	

	graph = TH1F("energyjet", "energyjet", 100, 0., 200.)
	graph2 = TH1F("energycherjet", "energycherjet", 100, 0., 200.)
	graph3 = TH1F("energyscinjet", "energyscinjet", 100, 0., 200.)
	graphmass = TH1F("mass_jet", "mass_jet", 100, 0., 200.)

	graph4 = TH1F("energy", "energy", 100, 0., 200.)
	graph5 = TH1F("energycher", "energycher", 100, 0., 200.)
	graph6 = TH1F("energyscin", "energyscin", 100, 0., 200.)

	#loop over events
	for Event in range(tree.GetEntries()):	

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
		energycher = float(sum(Calib_BarrelR_VectorSignalsCher)+sum(Calib_BarrelL_VectorSignalsCher))
		threshold = 0.0 #(GeV)	

		if energy>70.:
			
			#event displays with signals (p.e.)
			#if Event < 1:
				#displayfile.cd()
				#ROOTHistograms.create_eventdisplay_scin("Jet", BarrelR_VectorSignals, BarrelL_VectorSignals, "signal"+str(Event)) 
				#ROOTHistograms.create_eventdisplay_cher("Jet", BarrelR_VectorSignalsCher, BarrelL_VectorSignalsCher, "signal"+str(Event))
		
			#event displays with energy (GeV)	
			if Event<10:
					displayfile.cd()
					ROOTHistograms.create_eventdisplay_scin("Jet_energy", Calib_BarrelR_VectorSignals, Calib_BarrelL_VectorSignals, "energy"+str(Event), threshold) 
					ROOTHistograms.create_eventdisplay_cher("Jet_energy", Calib_BarrelR_VectorSignalsCher, Calib_BarrelL_VectorSignalsCher, "energy"+str(Event), threshold)	
			
			inputparticles_scin = []
			inputparticles_cher = []

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
					inputparticles_cher.append(fastjet.PseudoJet(towercher.Px(), towercher.Py(), towercher.Pz(), towercher.E()))
			
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
					inputparticles_cher.append(fastjet.PseudoJet(towercher.Px(), towercher.Py(), towercher.Pz(), towercher.E()))

			jet_def = fastjet.JetDefinition(fastjet.ee_genkt_algorithm, 2*math.pi, 1.)
			
			clust_seq = fastjet.ClusterSequence(inputparticles_scin, jet_def)	

			print "Event: "+str(Event)+" energy (GeV): "+str(energy)+" n-jets: "+str(len(clust_seq.exclusive_jets(int(2))))+" truth: "+str(len(inputparticles_scin))
			
			clust_seq_cher = fastjet.ClusterSequence(inputparticles_cher, jet_def)

			jet1_scin = fastjet.sorted_by_E(clust_seq.exclusive_jets(int(2)))[0]
			jet2_scin = fastjet.sorted_by_E(clust_seq.exclusive_jets(int(2)))[1]

			jet1_cher = fastjet.sorted_by_E(clust_seq_cher.exclusive_jets(int(2)))[0]
			jet2_cher = fastjet.sorted_by_E(clust_seq_cher.exclusive_jets(int(2)))[1]

			print "DeltaR jet1_scin: "+str(jet1_scin.delta_R(jet1_cher))+" "+str(jet1_scin.delta_R(jet2_cher))

			c = 0.34 #chi factor

			jet1, jet2 = mergejet(jet1_scin, jet2_scin, jet1_cher, jet2_cher)

			graph.Fill(jet1.e()+jet2.e())
			graph3.Fill(jet1_scin.e()+jet2_scin.e())
			graph2.Fill(jet1_cher.e()+jet2_cher.e())
			j = jet1+jet2
			graphmass.Fill(j.m())
			
			graph4.Fill((energy-c*energycher)/(1.-c))
			graph5.Fill(energycher)
			graph6.Fill(energy)
	
	graph.Write()
	graph2.Write()
	graph3.Write()
	graph4.Write()
	graph5.Write()
	graph6.Write()
	graphmass.Write()
	
jetdisplay()