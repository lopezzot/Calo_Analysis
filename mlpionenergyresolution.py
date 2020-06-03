import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TCanvas, TPad, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
#import calibration
#import calibration2 as calibration
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Dense
from matplotlib import pyplot as plt
from kerastuner.tuners import RandomSearch
import kerastuner

def build_model():
	model = keras.Sequential([
		keras.layers.Dense(8, activation='relu', input_shape=((2,))),
		keras.layers.Dense(8, activation='relu'),
		keras.layers.Dense(8, activation='relu'),
		#keras.layers.Dense(10, activation='relu'),
		keras.layers.Dense(1)
		])
	'''
	lr_schedule = keras.optimizers.schedules.InverseTimeDecay(
		0.001,
		decay_steps=78552*1000,
		decay_rate=1,
		staircase=False
		)
	'''
	
	#optimizer = keras.optimizers.RMSprop()
	#optimizer2 = keras.optimizers.SGD()
	model.compile(loss='mae', optimizer=keras.optimizers.Adam(learning_rate=0.0001),metrics=["mae","mse"])
	return model
'''
#hypermodel
def build_model_hp(hp):
	model = keras.Sequential([
		keras.layers.Dense(hp.Int('input_units_0',
		min_value=2,
		max_value=10,
		step=2), activation='relu', input_shape=((2,))),
		])
	for i in range(hp.Int('n_layers', 1, 4)):  # adding variation of layers.
		model.add(Dense(hp.Int('input_units_'+str(i), min_value=2, max_value=10, step=2), activation="relu"))
	model.add(Dense(1))
	hp_learning_rate = hp.Choice('learning_rate', values = [1e-2,1e-3,1e-4,1e-5])
	#hp_optimizer = hp.Choice('optimizer', ["sgd", "adam", "RMSprop"])
	model.compile(loss='mae', optimizer=keras.optimizers.Adam(learning_rate = hp_learning_rate), metrics=["mae","mse"])
	return model
#end hypermodel
'''

def recenergy(name):
	outputfile = "MLPionEnergyRes"+str(name)
	#displayfile = TFile(outputfile+".root","RECREATE")
	'''
	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	Energy = array('d')
	energyfractionscin = array('d')
	energyfractioncher = array('d')
	energyfraction = array('d')
	resolutionscin = array('d')
	resolutioncher = array('d')
	resolution = array('d')
	'''

	trainvector = np.array([[0,0]])
	labelvector = np.array([[0]])
	
	inputfiles = ["/home/lorenzo/Calo/results/trainPion_6_4_2020/trainPionFTFPBERT.root"]
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		#loop over events
		for Event in range(200000):	

			tree.GetEntry(Event)	
			if Event%10000==0:
				print "-> for training: ",Event 

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


			if PrimaryParticleEnergy/1000.<12.:
				cutleak = 0.5
			if PrimaryParticleEnergy/1000.>12. and PrimaryParticleEnergy/1000.<50.:
				cutleak = 1.0
			if PrimaryParticleEnergy/1000.>50.:
				cutleak = 3.0
			
			if (Leak/1000.+NeutrinoLeak/1000.)<cutleak and PrimaryParticleEnergy/1000. > 3.0:
				signalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
				signalcher = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)	
				#e_c = float(t[counter])#-(Leak/1000.+NeutrinoLeak/1000.)
				e_c = PrimaryParticleEnergy/1000.
				vector = np.array([[signalscin,signalcher]])
				label = np.array([[e_c]])
				trainvector = np.append(trainvector, vector, axis = 0)
				labelvector = np.append(labelvector, label, axis = 0)

	trainvector = np.delete(trainvector, 0, 0) #cancel first entry because it was set to 0,0,0
	labelvector = np.delete(labelvector, 0, 0)
	indices = np.arange(trainvector.shape[0])
	np.random.shuffle(indices)
	trainvector = trainvector[indices]
	labelvector = labelvector[indices] #random shuffle entries
	#normalize entries to maximum
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
		labelvector[i][0] = labelvector[i][0]/200. #normalize energy to 200. i.e. max energy

	trainvector2 = trainvector[:len(trainvector)/2]
	labelvector2 = labelvector[:len(labelvector)/2]
	evalvector = trainvector[len(trainvector)/2:len(trainvector)]
	evallabelvector = labelvector[len(labelvector)/2:len(labelvector)]
	print len(trainvector2), len(evalvector)
	return trainvector2, labelvector2, evalvector, evallabelvector, maxs, maxc

def gettestdata(name, t1, maxs, maxc):
	testvector = np.array([[0,0]])
	testlabelvector = np.array([[0]])
	
	inputfiles1 = ["/home/lorenzo/Calo/results/Pion_25_3_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t1]
	for counter, inputfile in enumerate(inputfiles1):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		#loop over events
		for Event in range(0,50000):	

			tree.GetEntry(Event)	
			if Event%10000==0:
				print "-> for predicting: ",Event 

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
			
			if float(t1[counter])<12.:
				cutleak = 0.5
			if float(t1[counter])>12. and float(t1[counter])<=50.:
				cutleak = 1.0
			if float(t1[counter])>50.:
				cutleak = 3.0

			if (Leak/1000.+NeutrinoLeak/1000.)<cutleak:
				signalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
				signalcher = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)	
				e_c = float(t1[counter])#-(Leak/1000.+NeutrinoLeak/1000.)
				vector = np.array([[signalscin,signalcher]])
				label = np.array([[e_c]])
				testvector = np.append(testvector, vector, axis = 0)
				testlabelvector = np.append(testlabelvector, label, axis = 0)

	testvector = np.delete(testvector, 0, 0)
	testlabelvector = np.delete(testlabelvector, 0, 0)
	indices = np.arange(testvector.shape[0])
	np.random.shuffle(indices)
	testvector = testvector[indices]
	testlabelvector = testlabelvector[indices]
	for i in range(len(testvector)):
		testvector[i][0] = testvector[i][0]/maxs
		testvector[i][1] = testvector[i][1]/maxc
		testlabelvector[i][0] = testlabelvector[i][0]/200.
	
	return testvector, testlabelvector

def doml(nepochs,trainvector, labelvector, evalvector, evallabelvector):
	'''
	#for tuning
	tuner = RandomSearch(build_model_hp,
						 objective = 'mse', 
						max_trials = 100,
						executions_per_trial = 3,
						directory = 'my_dir')
	
	tuner.search_space_summary()
	tuner.search(trainvector, labelvector, verbose = 0, epochs = 20, validation_data = (evalvector, evallabelvector))
	print len(trainvector), len(evalvector)
	#best_hps = tuner.get_best_hyperparameters(num_trials = 1)[0]
	#print best_hps
	print(tuner.results_summary())
	tuner.get_best_models()[0].save("hpmodel")
	#end tuning
	'''
	model = build_model()
	model.summary()
	
	history = model.fit(trainvector, labelvector, epochs = nepochs, validation_data = (evalvector,evallabelvector))
	
	history_dict = history.history
	train_mse = history_dict['mse']
	train_mse_array = array('d')
	for r in train_mse:
		train_mse_array.append(r)
	ep = range(1, len(train_mse)+1)
	epochs_array = array('d')
	for r in ep:
		epochs_array.append(r)
	eval_mse = history_dict['val_mse']
	eval_mse_array = array('d')
	for r in eval_mse:
		eval_mse_array.append(r)
	train_mse_graph = TGraph(len(epochs_array), epochs_array, train_mse_array)
	train_mse_graph.SetName("TrainMSE")
	eval_mse_graph = TGraph(len(epochs_array), epochs_array, eval_mse_array)
	eval_mse_graph.SetName("EvalMSE")

	#scatter plot
	scatterplot_eval = TH2F("scatteplot_eval", "scatteplot_eval", 400, 0.0, 250.0, 400, 0.0, 250.0)
	result_eval = model.predict(evalvector)
	for x in range(len(result_eval)):
		scatterplot_eval.Fill(evallabelvector[x]*200, result_eval[x]*200.)
	
	outputfile = "MLPionEnergyRes"+str(name)
	displayfile = TFile(outputfile+".root", "UPDATE")
	displayfile.cd()
	scatterplot_eval.Write()
	eval_mse_graph.SetMinimum(0.0001)
	eval_mse_graph.SetMaximum(0.001)
	train_mse_graph.SetMinimum(0.0001)
	train_mse_graph.SetMaximum(0.001)
	train_mse_graph.Write()
	eval_mse_graph.Write()

	return model

def mlpredict(model, i, testvector, testlabelvector):
	out = model.predict(testvector)
	scatplot = TH2F("scatplot_"+str(i), "scatplot_"+str(i), 400, 0.0, 200., 400, 0.0, 200.)
	hist = TH1F("hist"+str(i), "hist"+str(i), 400, 0.0, 200.)
	for i in range(len(out)):
		scatplot.Fill(testlabelvector[i]*200., out[i]*200.)
		#print testlabelvector[i]*200., out[i]*200.
		hist.Fill(out[i]*200.)
	
	outputfile = "MLPionEnergyRes"+str(name)
	displayfile = TFile(outputfile+".root", "UPDATE")
	displayfile.cd()
	scatplot.Write()
	
	hist.Fit("gaus")
	hist.Write()
	mean = hist.GetFunction("gaus").GetParameter(1)
	sigma = hist.GetFunction("gaus").GetParameter(2)
	sigmaerror = hist.GetFunction("gaus").GetParError(2)
	number_entries = hist.GetEntries()
		
	return mean, sigma, number_entries, sigmaerror

#names = ["FTFPBERT"]#, "FTFPBERT", "QGSPBERT", "QBBC"]
#names = ["FTFPBERTTRV", "QGSPBERT", "QBBC"]
#name = raw_input("physics list: ")
names = ["FTFPBERT"]

for name in names:
	trainvector, labelvector, evalvector, evallabelvector, maxs, maxc = recenergy(name)
	#model = doml(400,trainvector, labelvector, evalvector, evallabelvector)
	#model.save("saved_models/my_model")
	model = tf.keras.models.load_model("model_5")
	energy = [10,30,50,70,100,120,140,150]
	sqrtenergies = array('d')
	trueenergies = array('d')
	for e in energy:
		sqrtenergies.append(1./e**0.5)
		trueenergies.append(e)

	mean = array('d')
	meanerror = array('d')
	sigma_energy = array('d')
	sigma_energy_error = array('d')
	ratio = array('d')
	ratioerror = array('d')
	zeros = array('d')
	
	for e in energy:
		testvector, testlabelvector = gettestdata(name, [str(e)], maxs, maxc)
		m, s, number_entries , sigma_error =	mlpredict(model, e, testvector, testlabelvector)
		mean.append(m)
		meanerror.append(3*s/(float(number_entries)**0.5))
		sigma_energy.append(s/m)
		zeros.append(0.0)
		sigma_energy_e = ((s/(float(number_entries)**0.5))/m + sigma_error/s)*(s/m)
		sigma_energy_error.append(sigma_energy_e)

	for counter, e in enumerate(energy):
		ratio.append(mean[counter]/float(e))
		ratioerror.append(3*(s/(float(number_entries)**0.5)) / float(e))

	resgraph = TGraphErrors(len(sigma_energy), sqrtenergies, sigma_energy, zeros, sigma_energy_error)
	resgraph.SetName("resgraph")
	lingraph = TGraphErrors(len(trueenergies), trueenergies, mean, zeros, meanerror)
	lingraph.SetName("lingraph")
	ratiograph = TGraphErrors(len(trueenergies), trueenergies, ratio, zeros, ratioerror)
	ratiograph.SetName("ratiograph")
	displayfile = TFile("MLPionEnergyResFTFPBERT.root", "UPDATE")
	displayfile.cd()
	resgraph.Write()
	lingraph.Write()
	ratiograph.Write()

	ratiograph.SetMinimum(0.976)
	ratiograph.SetMaximum(1.024)
	ratiograph.GetXaxis().SetLimits(0.0,155.0)
	ratiograph.SetTitle("")
	ratiograph.GetXaxis().SetLabelSize(.105)
	ratiograph.GetYaxis().SetLabelSize(0.105)
	ratiograph.GetXaxis().SetNdivisions(520)
	ratiograph.GetYaxis().SetNdivisions(504)
	lingraph.SetTitle("")
	lingraph.SetMinimum(0.0)
	lingraph.SetMaximum(155.0)
	lingraph.GetXaxis().SetLimits(0.0,155.0)
	lingraph.GetXaxis().SetLabelSize(0.)
	lingraph.GetXaxis().SetNdivisions(520)
	lingraph.GetYaxis().SetNdivisions(520)
	c = TCanvas("c", "canvas", 800, 800)
	# Upper histogram plot is pad1
	pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
	pad1.SetBottomMargin(0.02)  # joins upper and lower plot
	pad1.SetLeftMargin(0.1)
	pad1.SetRightMargin(0.1)
	pad1.Draw()
	# Lower ratio plot is pad2
	c.cd()  # returns to main canvas before defining pad2
	pad2 = TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
	pad2.SetTopMargin(0.05)  # joins upper and lower plot
	pad2.SetBottomMargin(0.3)
	pad2.SetLeftMargin(0.1)
	pad2.SetRightMargin(0.1)
	pad2.Draw()
	pad1.cd()
	lingraph.Draw()
	pad2.cd()
	ratiograph.Draw("AP")
	ratioline = TF1("ratioline", str(np.mean(ratio)), 0., 160.)
	ratioline.SetLineColor(1)
	ratioline.SetLineWidth(1)
	ratioline.SetLineStyle(9)
	ratioline.Draw("same")
	ratioline.Write()
	c.Update()
	#c.SaveAs("MLratio.pdf")
	c.Write()
	