import map
import ROOTHistograms
from ROOT import TPad,gPad,TGaxis, TCanvas, TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
#import calibration
import calibration2 as calibration

machine = raw_input("On what machine running? (mac, linux, office) ")
if machine == "mac":
		path = str("/Users/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/Users/lorenzo/cernbox/work/Git-to-Mac/IDEA_Calorimeter_Union_data/")
if machine == "linux":
		path = str("/home/lorenzo/cernbox/work/Git-to-Mac/AnalysisIDEACalorimeter/")
		datapath = str("/home/lorenzo/Desktop/Calo/results/newenergyscan3/")
if machine == "office":
		datapath = str("/home/software/Calo/results/newenergyscan3_noangsmearing/")
'''
def reversexaxis(graph):
	graph.GetXaxis().SetLabelOffset(999)
	graph.GetXaxis().SetTickLength(0)
	gPad.Update()	
	newaxis = TGaxis(gPad.GetUxmax(),
                     gPad.GetUymin(),
                     gPad.GetUxmin(),
                     gPad.GetUymin(),
                     graph.GetXaxis().GetXmin(),
                     graph.GetXaxis().GetXmax(),
                     510,"-")
	newaxis.SetLabelOffset(-0.03)
	newaxis.Draw()
'''
def recenergy():
	outputfile = "EMEnergyRes"
	displayfile = TFile(outputfile+".root","RECREATE")

	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	Energy = array('d')
	energyfractionscin = array('d')
	ratio = array('d')
	ratioerror = array('d')
	energyfractionscinerror = array('d')
	energyfractioncher = array('d')
	energyfractionchererror = array('d')
	energyfraction = array('d')
	energyfractionerror = array('d')
	zeros = array('d')
	resolutionscin = array('d')
	resolutioncher = array('d')
	resolution = array('d')
	resolutionscinerror = array('d')
	resolutionchererror = array('d')
	resolutionerror = array('d')
	scin_sqrtenergies = array('d')
	cher_sqrtenergies = array('d')
	
	#inputfiles = sorted(glob.glob(datapath+"*"), key=os.path.getmtime) #get files from tower 1 to 75 ordered by creation time
	t = [5, 10,40, 60, 80, 100, 150]
	t = [10,40,60,80,100]
	t = [10,30,50,70,90,100,110,140,150,250]
	energies = array('d', t)
	sqrtenergies = array('d',[1/(x**0.5) for x in energies])
	inputfiles = ["/home/software/Calo/results/newenergyscan4/Electron_"+str(i)+".root" for i in t]
	inputfiles = ["/home/software/Calo/results/newresults/Electron_11_4_2020/Electron_"+str(i)+".root" for i in t]	
	#newenergyscan3_noangsmearing per linearita da 5 a 150 gev
	#newenergyscan4 per risoluzione
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		if t[counter]<35.:
			ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 200, 0., 40.)
			CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 200, 0., 40.)	
			RecEnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 200, 0., 40.)
		if t[counter]>35. and t[counter]<400.:
			ScinEnergyHist = TH1F("scinenergy_", str(counter+1)+"_scin", 1500, 0., 300.)
			CherEnergyHist = TH1F("cherenergy_", str(counter+1)+"_cher", 1500, 0., 300.)	
			RecEnergyHist = TH1F("RecEnergy_",str(counter+1)+"_Energy", 1500, 0., 300.)

		Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 3000, 0., 30000)
		EnergyHist = TH1F("Energy_",str(counter+1)+"_Energy", 500, 0., 200.)
		
		scatterplot = TH2F("scatterplot_", str(counter+1), int(800), 0., 200., int(800), 0., 200.)

		#loop over events
		for Event in range(10000):	

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
			#sigmascin = 0.15*(energyscin**0.5)+0.012*energyscin #old one
			sigmascin = 0.177*(energyscin**0.5)+0.006*energyscin
			CherEnergyHist.Fill(energycher)
			#sigmacher = 0.18*(energycher**0.5)+0.0045*energycher #old one
			sigmacher = 0.194*(energycher**0.5)+0.001*energycher
			RecEnergyHist.Fill((energyscin/(sigmascin**2)+energycher/(sigmacher**2))/(1/sigmascin**2+1/sigmacher**2))

			scatterplot.Fill(energyscin, energycher)
		
		print energies[counter], ScinEnergyHist.GetMean(), CherEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		ScinEnergyHist.Fit("gaus")
		CherEnergyHist.Fit("gaus")
		RecEnergyHist.Fit("gaus")
		RecEnergyHist.Write()
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		#Signalscinhist.Write()
		#EnergyHist.Write()
		scatterplot.Write()
		scin_sqrtenergies.append(1./(ScinEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		cher_sqrtenergies.append(1./(CherEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		MeanEnergyScin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		MeanEnergyCher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1))
		
		energyfractionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		energyfractionscinerror.append(3.0*ScinEnergyHist.GetFunction("gaus").GetParError(1))		
		energyfractioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(1))
		energyfractionchererror.append(3.0*CherEnergyHist.GetFunction("gaus").GetParError(1))
		energyfraction.append(RecEnergyHist.GetFunction("gaus").GetParameter(1))
		energyfractionerror.append(3.0*RecEnergyHist.GetFunction("gaus").GetParError(1))
		zeros.append(0.0)		
		ratio.append(RecEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		z = 3.0*RecEnergyHist.GetFunction("gaus").GetParError(1)/RecEnergyHist.GetFunction("gaus").GetParameter(1)
		#z = RecEnergyHist.GetFunction("gaus").GetParameter(2)/100./RecEnergyHist.GetFunction("gaus").GetParameter(1)
		ratioerror.append(z*RecEnergyHist.GetFunction("gaus").GetParameter(1)/energies[counter])
		
		resolutionscin.append(ScinEnergyHist.GetFunction("gaus").GetParameter(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		z = ScinEnergyHist.GetFunction("gaus").GetParError(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(2)+ScinEnergyHist.GetFunction("gaus").GetParError(1)/ScinEnergyHist.GetFunction("gaus").GetParameter(1)
		resolutionscinerror.append(z*ScinEnergyHist.GetFunction("gaus").GetParameter(2)/ScinEnergyHist.GetFunction("gaus").GetParameter(1))
		resolutioncher.append(CherEnergyHist.GetFunction("gaus").GetParameter(2)/CherEnergyHist.GetFunction("gaus").GetParameter(1))
		z = CherEnergyHist.GetFunction("gaus").GetParError(2)/CherEnergyHist.GetFunction("gaus").GetParameter(2)+CherEnergyHist.GetFunction("gaus").GetParError(1)/CherEnergyHist.GetFunction("gaus").GetParameter(1)
		resolutionchererror.append(z*CherEnergyHist.GetFunction("gaus").GetParameter(2)/CherEnergyHist.GetFunction("gaus").GetParameter(1))		
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))		
		z = RecEnergyHist.GetFunction("gaus").GetParError(2)/RecEnergyHist.GetFunction("gaus").GetParameter(2)+RecEnergyHist.GetFunction("gaus").GetParError(1)/RecEnergyHist.GetFunction("gaus").GetParameter(1)		
		resolutionerror.append(z*RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))		
		
	ratioGraph = TGraphErrors(len(energies), energies, ratio, zeros, ratioerror)
	ratioGraph.SetName("ratio")
	ratioGraph.SetMinimum(0.976)
	ratioGraph.SetMaximum(1.024)
	ratioGraph.GetXaxis().SetLimits(0.0,155.0)
	ratioGraph.SetTitle("")
	ratioGraph.GetXaxis().SetLabelSize(.105)
	ratioGraph.GetYaxis().SetLabelSize(0.105)
	ratioGraph.GetXaxis().SetNdivisions(520)
	ratioGraph.GetYaxis().SetNdivisions(504)
	ratioGraph.Write()
	LinearityGraph = TGraphErrors(len(energies), energies, energyfraction, zeros, energyfractionerror)
	LinearityGraph.SetName("LinearityGraph")
	LinearityGraph.Write()
	LinearityGraph.SetTitle("")
	LinearityGraph.SetMinimum(0.0)
	LinearityGraph.SetMaximum(155.0)
	LinearityGraph.GetXaxis().SetLimits(0.0,155.0)
	LinearityGraph.GetXaxis().SetLabelSize(0.)
	LinearityGraph.GetXaxis().SetNdivisions(520)
	LinearityGraph.GetYaxis().SetNdivisions(520)
	LinearityGraphScin = TGraphErrors(len(energies), energies, energyfractionscin, zeros, energyfractionscinerror)
	LinearityGraphCher = TGraphErrors(len(energies), energies, energyfractioncher, zeros, energyfractionchererror)
	LinearityGraphCher.SetName("LinearityGraphCher")
	LinearityGraphCher.Write()
	LinearityGraphScin.SetName("LinearityGraphScin")
	LinearityGraphScin.Write()

	ResolutionGraphScin = TGraphErrors(len(energies), scin_sqrtenergies, resolutionscin, zeros, resolutionscinerror)
	func = TF1("func", "[0]*x+[1]", 0.05, 0.40)
	ResolutionGraphCher = TGraphErrors(len(energies), cher_sqrtenergies, resolutioncher, zeros, resolutionchererror)
	ResolutionGraphScin.Fit("func", "R")
	ResolutionGraphCher.Fit("func", "R")
	ResolutionGraphScin.SetName("ResolutionGraphScin")
	ResolutionGraphScin.GetXaxis().SetLimits(0.0,0.6)
	ResolutionGraphScin.SetMinimum(0.0)
	ResolutionGraphScin.SetMaximum(0.1)
	ResolutionGraphScin.Write()
	ResolutionGraphCher.SetName("ResolutionGraphCher")
	ResolutionGraphCher.GetXaxis().SetLimits(0.0,0.6)
	ResolutionGraphCher.SetMinimum(0.0)
	ResolutionGraphCher.SetMaximum(0.1)
	ResolutionGraphCher.Write()
	ResolutionGraph = TGraphErrors(len(energies), sqrtenergies, resolution, zeros, resolutionerror)
	ResolutionGraph.Fit("func", "R")
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.GetXaxis().SetLimits(0.0,0.6)
	ResolutionGraph.SetMinimum(0.0)
	ResolutionGraph.SetMaximum(0.1)
	ResolutionGraph.Write()

	'''
	rd52copper = array('d', [0.04478505426185217, 0.027392527130926082, 0.02420093893609386, 0.02229837387624884, 0.020999999999999998])
	rd52graph = TGraph(len(energies), sqrtenergies, rd52copper)
	rd52graph.SetName("rd52resolution")
	rd52graph.Write()


	EMResolutions = TMultiGraph()
	EMResolutions.Add(ResolutionGraphScin)
	EMResolutions.Add(ResolutionGraphCher)
	EMResolutions.Add(ResolutionGraph)
	EMResolutions.Add(rd52graph)
	EMResolutions.SetName("EMResolutions")
	EMResolutions.Write()
	'''
	#Linearities = TMultiGraph()
	
	#Linearities.Add(LinearityGraph)
	#Linearities.Add(LinearityGraphScin)
	#Linearities.Add(LinearityGraphCher)
	#Linearities.SetName("Linearities")
	#Linearities.Write()
	
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
	LinearityGraph.Draw()
	pad2.cd()
	ratioGraph.Draw()
	c.Update()
	c.SaveAs("test.pdf")
	c.Write()


	ca = TCanvas("ca", "canvas1", 800, 800)
	ca.cd()
	gPad.DrawFrame(0,0,0.6,0.1)
	fa=TF1("fa","exp(x)",0,2)
	Axis=TGaxis(0,0.1,0.6,0.1,"fa",510,"-")
	Axis.SetFunction("fa")
	Axis.Draw()
	#ResolutionGraph = TGraphErrors(len(energies), sqrtenergies, resolution, zeros, resolutionerror)
	#ResolutionGraph.Fit("func", "R")
	#ResolutionGraph.SetName("ResolutionGraph")
	#ResolutionGraph.GetXaxis().SetLimits(0.0,0.6)
	#ResolutionGraph.SetMinimum(0.0)
	#ResolutionGraph.SetMaximum(0.1)

	#ResolutionGraph.Draw()
	ca.Write()
	
	


recenergy()