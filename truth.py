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
import newmap_truth
	
def mergejet(jet1_scin, jet2_scin, jet1_cher, jet2_cher):
	deltar1 = jet1_scin.delta_R(jet1_cher)
	deltar2 = jet1_scin.delta_R(jet2_cher)

	jet1 = None
	jet2 = None

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

#used to get w mass by subtracting the muon from truth
def jetdisplay():

	inputfile1 = "wwlj_truth.root"
	inputfile2 = "wwlj.root"

	inputfile1 = TFile(inputfile1)
	inputfile2 = TFile(inputfile2)
	print "Analyzing: "+str(inputfile1)+" \n"
	tree1 = TTree()
	tree2 = TTree()
	inputfile1.GetObject("truth", tree1)	
	inputfile2.GetObject("B4", tree2)	
	tree1.AddFriend(tree2)

	outputfile = "wwlj_output"
	displayfile = TFile(outputfile+".root","RECREATE")

	graphmass = TH1F("mass_jet", "mass_jet", 100, 0., 200.)
	graphmass_truth = TH1F("mass_jet_truth", "mass_jet_truth", 100, 0., 200.)
	
	#loop over events
	for Event in range(int(10)):	

		tree1.GetEntry(Event)	

		#Set values of the tree
		numtru=tree1.mcs_n
		print numtru

		muvec = []
		inputparticles_tru = []
		nmuon=0
		#loop over true particles
		for itru in range(0,numtru):
			partid = tree1.mcs_pdgId[itru]
				   
			#for particle depositing in calo, store them as input for jet building
			if abs(partid) != 13 and  abs(partid) !=12 and abs(partid) != 14 and abs(partid) != 16:
				trup= TLorentzVector()
				trup.SetPtEtaPhiM(tree1.mcs_pt[itru], tree1.mcs_eta[itru], tree1.mcs_phi[itru], tree1.mcs_m[itru])
				inputparticles_tru.append(fastjet.PseudoJet(trup.Px(), trup.Py(), trup.Pz(), trup.E()))
			#store muons in event
			if abs(partid)==13: 
				muon=TLorentzVector()
				muon.SetPtEtaPhiM(tree1.mcs_pt[itru], tree1.mcs_eta[itru], tree1.mcs_phi[itru], tree1.mcs_m[itru])
				muvec.append(muon)
				nmuon=nmuon+1
		print " nmuon ",nmuon

		#now build truth jets
		jet_def = fastjet.JetDefinition(fastjet.ee_genkt_algorithm, 2*math.pi, 1.)
		clust_seq = fastjet.ClusterSequence(inputparticles_tru, jet_def)
		jetexc = fastjet.sorted_by_E(clust_seq.exclusive_jets(int(2)))
		print "*********** jets ************"
		for jet in jetexc: 
			print jet.e(), jet.eta(), jet.phi()
		print "*********** muons ************"
		for muon in muvec:
			print muon.E(), muon.Eta(), muon.Phi()
		
		jet1_truth = jetexc[0]
		jet2_truth = jetexc[1]
		j = jet1_truth+jet2_truth
		graphmass_truth.Fill(j.m())

		# now handle calo sim
		BarrelR_VectorSignals = tree2.VectorSignalsR
		BarrelL_VectorSignals = tree2.VectorSignalsL 
		BarrelR_VectorSignalsCher = tree2.VectorSignalsCherR
		BarrelL_VectorSignalsCher = tree2.VectorSignalsCherL
		VectorR = tree2.VectorR
		VectorL = tree2.VectorL

		Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
		Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
		Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
		Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)

		energy = float(sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals))
		print " simulated energy ", energy
		if(energy>0):
			threshold=0.1
			inputparticles_scin = []
			inputparticles_cher = []
				   
			#right part
			for towerindex in range(75*36):
				theta, phi, eta = newmap_truth.maptower(towerindex, "right")
				energy_scin = Calib_BarrelR_VectorSignals[towerindex]
				pt_scin = energy_scin*np.sin(theta*math.pi/180.)

				energy_cher = Calib_BarrelR_VectorSignalsCher[towerindex]
				pt_cher = energy_cher*np.sin(theta*math.pi/180.)

				towerscin = TLorentzVector()
				towerscin.SetPtEtaPhiM(pt_scin, eta, phi*math.pi/180., 0.)
				towercher = TLorentzVector()
				towercher.SetPtEtaPhiM(pt_cher, eta, phi*math.pi/180., 0.)
				deltamumin=999999.
				for muon in muvec:
					deltaR=abs(towerscin.DeltaR(muon))
					if deltaR<deltamumin:
						deltamumin=deltaR 
				if energy_scin > threshold:
					if deltamumin<0.1:
						print " deltamumin ", deltamumin
					if deltamumin>0.1:
						inputparticles_scin.append(fastjet.PseudoJet(towerscin.Px(), towerscin.Py(), towerscin.Pz(), towerscin.E()))
						inputparticles_cher.append(fastjet.PseudoJet(towercher.Px(), towercher.Py(), towercher.Pz(), towercher.E()))
			
			#left part	
			for towerindex in range(75*36):
				theta, phi, eta = newmap_truth.maptower(towerindex, "left")
				energy_scin = Calib_BarrelL_VectorSignals[towerindex]
				pt_scin = energy_scin*np.sin(theta*math.pi/180.)

				energy_cher = Calib_BarrelL_VectorSignalsCher[towerindex]
				pt_cher = energy_cher*np.sin(theta*math.pi/180.)

				towerscin = TLorentzVector()
				towerscin.SetPtEtaPhiM(pt_scin, eta, phi*math.pi/180., 0.)
				towercher = TLorentzVector()
				towercher.SetPtEtaPhiM(pt_cher, eta, phi*math.pi/180., 0.)
				deltamumin=999999.
				for muon in muvec:
					deltaR=abs(towerscin.DeltaR(muon))
					if deltaR<deltamumin:
						deltamumin=deltaR 
				if energy_scin > threshold:
					if deltamumin<0.1:
						print " deltamumin ", deltamumin
					if deltamumin>0.1:
						inputparticles_scin.append(fastjet.PseudoJet(towerscin.Px(), towerscin.Py(), towerscin.Pz(), towerscin.E()))
						inputparticles_cher.append(fastjet.PseudoJet(towercher.Px(), towercher.Py(), towercher.Pz(), towercher.E()))
		
		print "len: ",len(inputparticles_scin)
		print "lencher: ",len(inputparticles_cher)

		jet_def = fastjet.JetDefinition(fastjet.ee_genkt_algorithm, 2*math.pi, 1.)
		
		clust_seq = fastjet.ClusterSequence(inputparticles_scin, jet_def)	
		clust_seq_cher = fastjet.ClusterSequence(inputparticles_cher, jet_def)

		print "n jet: ", len(clust_seq.exclusive_jets(int(2))), len(clust_seq_cher.exclusive_jets(int(2)))
		jet1_scin = clust_seq.exclusive_jets(int(2))[0]
		jet2_scin = clust_seq.exclusive_jets(int(2))[1]

		jet1_cher = clust_seq_cher.exclusive_jets(int(2))[0]
		jet2_cher = clust_seq_cher.exclusive_jets(int(2))[1]

		#merge jet
		jet1, jet2 = mergejet(jet1_scin, jet2_scin, jet1_cher, jet2_cher)
		jet = jet1+jet2
		graphmass.Fill(jet.m())
	graphmass.Write()
	graphmass_truth.Write()


jetdisplay()
