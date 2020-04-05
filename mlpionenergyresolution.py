import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
#import calibration
#import calibration2 as calibration
import tensorflow as tf
from tensorflow import keras
from matplotlib import pyplot as plt

def build_model():
	model = keras.Sequential([
		keras.layers.Dense(40, activation='relu', input_shape=((2,))),
		keras.layers.Dense(40, activation='relu'),
		keras.layers.Dense(40, activation='relu'),
		keras.layers.Dense(1)
		])
	optimizer = keras.optimizers.RMSprop(0.001)
	model.compile(loss='mse', optimizer=optimizer,metrics=['mae', "mse"])
	return model

def recenergy(name):
	outputfile = "MLPionEnergyRes"+str(name)
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

	trainvector = np.array([[0,0]])
	labelvector = np.array([[0]])
	testvector = np.array([[0,0]])
	testlabelvector = np.array([[0]])
	
	t = [10,30,50,70,100,120,140,150]
	energies = array('d',t)
	inputfiles = ["/home/lorenzo/Calo/results/Pion_25_3_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t]
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_"+str(t[counter]), str(t[counter])+"_scin", 400, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_"+str(t[counter]), str(t[counter])+"_cher", 400, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_"+str(t[counter]),str(t[counter])+"_Energy", 400, 0., 200.)
		
		#loop over events
		for Event in range(10000):	

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

			if (Leak/1000.+NeutrinoLeak/1000.)<3.0:
				signalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
				signalcher = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)	
				e_c = float(t[counter])#-(Leak/1000.+NeutrinoLeak/1000.)
				vector = np.array([[signalscin,signalcher]])
				label = np.array([[e_c]])
				trainvector = np.append(trainvector, vector, axis = 0)
				labelvector = np.append(labelvector, label, axis = 0)

	t1 = [10, 30, 50, 70, 100, 120, 140, 150]
	t1 = [100]
	energies = array('d',t1)
	inputfiles1 = ["/home/lorenzo/Calo/results/Pion_25_3_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t1]
	for counter, inputfile in enumerate(inputfiles1):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		#loop over events
		for Event in range(10000,20000):	

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

			if (Leak/1000.+NeutrinoLeak/1000.)<3.0:
				signalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
				signalcher = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)	
				e_c = float(t1[counter])#-(Leak/1000.+NeutrinoLeak/1000.)
				vector = np.array([[signalscin,signalcher]])
				label = np.array([[e_c]])
				testvector = np.append(testvector, vector, axis = 0)
				testlabelvector = np.append(testlabelvector, label, axis = 0)

	trainvector = np.delete(trainvector, 0, 0)
	labelvector = np.delete(labelvector, 0, 0)
	indices = np.arange(trainvector.shape[0])
	np.random.shuffle(indices)
	trainvector = trainvector[indices]
	labelvector = labelvector[indices]
	maxs = 0
	maxc = 0
	max1 = 0
	max2 = 0
	for i in range(len(trainvector)):
		max1 = trainvector[i][0]
		if max1 > maxs:
			maxs = max1
		max2 = trainvector[i][1]
		if max2 > maxc:
			maxc = max2
	for i in range(len(trainvector)):
		trainvector[i][0] = trainvector[i][0]/maxs
		trainvector[i][1] = trainvector[i][1]/maxc
		labelvector[i][0] = labelvector[i][0]/150.

	testvector = np.delete(testvector, 0, 0)
	testlabelvector = np.delete(testlabelvector, 0, 0)
	indices = np.arange(testvector.shape[0])
	np.random.shuffle(indices)
	testvector = testvector[indices]
	testlabelvector = testlabelvector[indices]
	for i in range(len(testvector)):
		testvector[i][0] = testvector[i][0]/maxs
		testvector[i][1] = testvector[i][1]/maxc
		testlabelvector[i][0] = testlabelvector[i][0]/150.

	return trainvector, labelvector, testvector, testlabelvector

#names = ["FTFPBERT"]#, "FTFPBERT", "QGSPBERT", "QBBC"]
#names = ["FTFPBERTTRV", "QGSPBERT", "QBBC"]
#name = raw_input("physics list: ")
names = []
names.append("FTFPBERT")

for name in names:
	trainvector, labelvector, testvector, testlabelvector = recenergy(name)
	model = build_model()
	model.summary()
	history = model.fit(trainvector, labelvector, epochs = 50)
	history_dict = history.history
	train_mse = history_dict['mse']
	ep = range(1, len(train_mse)+1)
	plt.plot(ep, train_mse, 'bo', label='Training_mse')
	plt.xlabel('epochs')
	plt.ylabel('MSE')
	plt.ylim([0,0.01])
	plt.show()
	
	out = model.predict(testvector).flatten()
	scatplot = TH2F("scatplot", "scatplot", 400, -0.5, 200., 400, -0.5, 200.)
	hist = TH1F("hist", "hist", 100, -10., 10.)
	for i in range(len(out)):
		scatplot.Fill(testlabelvector[i]*150., out[i]*150.)
		print testlabelvector[i]*150., out[i]*150.
		hist.Fill(out[i]*150.-testlabelvector[i]*150.)
	outputfile = "MLPionEnergyRes"+str(name)
	displayfile = TFile(outputfile+".root","RECREATE")
	displayfile.cd()
	scatplot.Write()
	hist.Write()
	print trainvector
	print testvector
	
	

	