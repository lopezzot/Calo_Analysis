import map
import ROOTHistograms
from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1
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
		datapath = str("/home/software/Calo/results/NewTowerScan2/")

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
	Energy = array('d')

	ResponseMeanScin = array('d')
	ResponseMeanCher = array('d')
	ResponseRMSScin = array('d')
	ResponseRMSCher = array('d')
	ResponseZeros = array('d')
	Responseerrorsscin = array('d')
	Responseerrorscher = array('d')
	energytot = array('d')
	scinsignaltot = array('d')
	chersignaltot = array('d')
	ScinSignalTot = array('d')
	

	inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	#inputfiles = ["/home/software/Calo/Calibrations/Calib_"+str(i)+".root" for i in range(1,76)]
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinHist = TH1F("scin_", str(counter+1)+"_scin", 1000, 0., 1000.)
		CherHist = TH1F("cher_", str(counter+1)+"_cher", 600, 0., 600.)	
		EnergyTowerHist = TH1F("Energy_",str(counter+1)+"_Energy", 200, 0., 200.)
		EnergyHist = TH1F("E_", str(counter+1)+"_Energy", 200, 0., 200.)
		ResponseScinHist = TH1F("responsescin_", str(counter+1)+"_scin", 1000, 0., 1000.)
		ResponseCherHist = TH1F("responsecher_", str(counter+1)+"_cher", 600, 0., 600.)

		Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 30000, 0., 30000)

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
			energytower = max(VectorR)/1000
			
			ScinHist.Fill(signalscin/energytower)	
			CherHist.Fill(signalcher/energytower)	
			EnergyTowerHist.Fill(energytower)

			totalsignalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
			totalsignalcher = sum(BarrelR_VectorSignalsCher)+sum(BarrelL_VectorSignalsCher)
			totalenergy = (sum(VectorR)+sum(VectorL))/1000
			EnergyHist.Fill(totalenergy)


			ResponseScinHist.Fill(totalsignalscin/totalenergy)
			ResponseCherHist.Fill(totalsignalcher/totalenergy)
			Signalscinhist.Fill(totalsignalscin)

			energytot.append(totalenergy)
			scinsignaltot.append(totalsignalscin)
			chersignaltot.append(totalsignalcher)

			if Event < 1:
				print "Max found at: "+str(list(BarrelR_VectorSignals).index(signalscin))+str(list(BarrelR_VectorSignalsCher).index(signalcher))+str(list(VectorR).index(energytower*1000))+" for file "+str(counter+1) #to check tower mostly hitten is the correct one
				displayfile.cd()
				ROOTHistograms.create_eventdisplay_scin(PrimaryParticleName, BarrelR_VectorSignals, BarrelL_VectorSignals, str(counter)) 
				ROOTHistograms.create_eventdisplay_cher(PrimaryParticleName, BarrelR_VectorSignalsCher, BarrelL_VectorSignalsCher, str(counter))

		print np.mean(energytot), np.mean(scinsignaltot), np.mean(chersignaltot)
		displayfile.cd()
		gStyle.SetOptStat(111)
		ScinHist.Fit("gaus")
		CherHist.Fit("gaus")
		ScinHist.Write()
		CherHist.Write()
		EnergyTowerHist.Write()
		MeanScin.append(ScinHist.GetFunction("gaus").GetParameter(1))
		MeanCher.append(CherHist.GetFunction("gaus").GetParameter(1))
		RMSScin.append(ScinHist.GetRMS())
		RMSCher.append(CherHist.GetRMS())
		EnergyTower.append(EnergyTowerHist.GetMean())
		Energy.append(EnergyHist.GetMean())
		ScinSignalTot.append(Signalscinhist.GetMean())
		Signalscinhist.Write()
		#print ScinHist.GetMean(), CherHist.GetMean(), EnergyTowerHist.GetMean()
		#print ScinHist.GetRMS(), CherHist.GetRMS()
		Tower.append(counter+1)
		Zeros.append(0.)
		errorsscin.append(ScinHist.GetRMS()/(3000**0.5))
		errorscher.append(CherHist.GetRMS()/(3000**0.5))

		ResponseScinHist.Fit("gaus")
		ResponseCherHist.Fit("gaus")
		ResponseScinHist.Write()
		ResponseCherHist.Write()
		ResponseMeanScin.append(ResponseScinHist.GetFunction("gaus").GetParameter(1))
		ResponseMeanCher.append(ResponseCherHist.GetFunction("gaus").GetParameter(1))
		ResponseRMSScin.append(ResponseScinHist.GetRMS())
		ResponseRMSCher.append(ResponseCherHist.GetRMS())
		Responseerrorsscin.append(ResponseScinHist.GetRMS()/(3000**0.5))
		Responseerrorscher.append(ResponseCherHist.GetRMS()/(3000**0.5))

	n = len(MeanCher)
	ScinTotGraph = TGraph(n, Tower, ScinSignalTot)
	ScinTotGraph.SetName("ScinTotGraph")
	ScinTotGraph.Write()
	EnergyGraph = TGraph(n, Tower, Energy)
	EnergyGraph.SetName("EnergyGraph")
	EnergyGraph.Write()
	EnergyTowerGraph = TGraph(n, Tower, EnergyTower)
	EnergyTowerGraph.SetName("EnergyTower")
	EnergyTowerGraph.Write()
	RMSGraphScin = TGraph(n, Tower, RMSScin)
	RMSGraphScin.SetName("Calibration_RMSGraphScin")
	RMSGraphCher = TGraph(n, Tower, RMSCher)
	RMSGraphCher.SetName("Calibration_RMSGraphCher")
	MeanGraphScin = TGraphErrors(n, Tower, MeanScin, Zeros, errorsscin)
	MeanGraphScin.SetName("Calibration_Scin")
	MeanGraphCher = TGraphErrors(n, Tower, MeanCher, Zeros, errorscher)
	MeanGraphCher.SetName("Calibration_Cher")
	x = array('d', (0., 90., 90., 0.))
	y = array('d', (np.mean(MeanScin)-0.01*np.mean(MeanScin), np.mean(MeanScin)-0.01*np.mean(MeanScin), np.mean(MeanScin)+0.01*np.mean(MeanScin), np.mean(MeanScin)+0.01*np.mean(MeanScin)))
	Fillgraph = TGraph(4, x, y )
	Fillgraph.SetName("Calibration_banscin")
	linefillgraph = TF1("CalibrationMeanScin", str(np.mean(MeanScin)), 0., 90.)
	linefillgraph.Write()
	Fillgraph.Write()
	x2 = array('d', (0., 90., 90., 0.))
	y2 = array('d', (np.mean(MeanCher)-0.01*np.mean(MeanCher), np.mean(MeanCher)-0.01*np.mean(MeanCher), np.mean(MeanCher)+0.01*np.mean(MeanCher), np.mean(MeanCher)+0.01*np.mean(MeanCher)))
	Fillgraph2 = TGraph(4, x2, y2 )
	Fillgraph2.SetName("Calibration_bancher")
	linefillgraph2 = TF1("CalibrationMeanCher", str(np.mean(MeanCher)), 0., 90.)
	linefillgraph2.Write()
	Fillgraph2.Write()
	MeanGraphCher.Write()
	MeanGraphScin.Write()
	RMSGraphCher.Write()
	RMSGraphScin.Write()

	ResponseRMSGraphScin = TGraph(n, Tower, ResponseRMSScin)
	ResponseRMSGraphScin.SetName("ResponseRMSGraphScin")
	ResponseRMSGraphCher = TGraph(n, Tower, ResponseRMSCher)
	ResponseRMSGraphCher.SetName("ResponseRMSGraphCher")
	ResponseMeanGraphScin = TGraphErrors(n, Tower, ResponseMeanScin, Zeros, Responseerrorsscin)
	ResponseMeanGraphScin.SetName("ResponseMeanGraphScin")
	ResponseMeanGraphCher = TGraphErrors(n, Tower, ResponseMeanCher, Zeros, Responseerrorscher)
	ResponseMeanGraphCher.SetName("ResponseMeanGraphCher")
	x = array('d', (0., 90., 90., 0.))
	y = array('d', (np.mean(ResponseMeanScin)-0.01*np.mean(ResponseMeanScin), np.mean(ResponseMeanScin)-0.01*np.mean(ResponseMeanScin), np.mean(ResponseMeanScin)+0.01*np.mean(ResponseMeanScin), np.mean(ResponseMeanScin)+0.01*np.mean(ResponseMeanScin)))
	Fillgraph = TGraph(4, x, y )
	Fillgraph.SetName("ResponseBan_scin")
	linefillgraph = TF1("ResponseMeanScin", str(np.mean(ResponseMeanScin)), 0., 90.)
	linefillgraph.Write()
	Fillgraph.Write()
	x2 = array('d', (0., 90., 90., 0.))
	y2 = array('d', (np.mean(ResponseMeanCher)-0.01*np.mean(ResponseMeanCher), np.mean(ResponseMeanCher)-0.01*np.mean(ResponseMeanCher), np.mean(ResponseMeanCher)+0.01*np.mean(ResponseMeanCher), np.mean(ResponseMeanCher)+0.01*np.mean(ResponseMeanCher)))
	Fillgraph2 = TGraph(4, x2, y2 )
	Fillgraph2.SetName("ResponseBan_cher")
	linefillgraph2 = TF1("ResponseMeanCher", str(np.mean(ResponseMeanCher)), 0., 90.)
	linefillgraph2.Write()
	Fillgraph2.Write()
	ResponseMeanGraphCher.Write()
	ResponseMeanGraphScin.Write()
	ResponseRMSGraphCher.Write()
	ResponseRMSGraphScin.Write()

towercalibration()