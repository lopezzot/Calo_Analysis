import map
import ROOTHistograms
from ROOT import TTree, TProfile, TFile, TCutG, TH1F, TCanvas, TPad, TGraph, TGraphErrors, gStyle, TH2F, TLine, TF1, TMultiGraph
import glob
from array import array
import os
import newmap
import numpy as np
#import calibration
import calibration2 as calibration

def recenergy(name, chivalue):
	outputfile = "KaonEnergyRes_"+str(chivalue)+"_"+str(name)
	displayfile = TFile(outputfile+".root","RECREATE")
	'''
	if name == "FTFPBERT":
		chi = 0.3
	if name == "FTFPBERTTRV":
		chi = 0.3
	if name == "QGSPBERT":
		chi = 0.3
	if name == "QBBC":
		chi = 0.3
	'''
	MeanEnergyScin = array('d')
	MeanEnergyCher = array('d')
	Energy = array('d')
	Energyerror = array('d')
	energyfractionscin = array('d')
	energyfractionscinerror = array('d')
	energyfractioncher = array('d')
	energyfractionchererror = array('d')
	energyfraction = array('d')
	energyfractionerror = array('d')
	resolutionscin = array('d')
	resolutionscinerror = array('d')
	resolutioncher = array('d')
	resolutionchererror = array('d')
	resolution = array('d')
	resolutionerror = array('d')
	chiarray = array('d')
	chierrorarray = array('d')
	zeros = array('d')
	containment = array('d')
	containmenterror = array('d')

	t = [10,30,50,70,100,120,140,150]
	t = [100]
	sqrtenergies = array('d',[1/(x**0.5) for x in t])
	energies = array('d',t)
	#inputfiles = ["/home/software/Calo/results/energycont_2p0m/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/software/Calo/results/pionenergyscan_QGSPBICHP/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Desktop/Calo/newresults/FTFPBERTTRV/Pion_"+str(i)+"_FTFPBERTTRV_office.root" for i in t]
	#inputfiles = ["/Users/lorenzo/Desktop/ToPC/newresults/"+str(name)+"/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Calo/results/Pion_25_3_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Calo/results/geant4.10.4.p01/Pion_30_4_2020/"+str(name)+""+"/Pion_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Calo/results/Proton_25_3_2020/"+str(name)+""+"/Proton_"+str(i)+".root" for i in t]
	#inputfiles = ["/home/lorenzo/Calo/results/Neutron_25_3_2020/"+str(name)+""+"/Neutron_"+str(i)+".root" for i in t]
	inputfiles = ["/home/lorenzo/Calo/results/Kaon_5_4_2020/"+str(name)+""+"/Kaon_"+str(i)+".root" for i in t]

	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("B4", tree)	

		ScinEnergyHist = TH1F("scinenergy_"+str(t[counter]), str(t[counter])+"_scin", 400, 0., 200.)
		CherEnergyHist = TH1F("cherenergy_"+str(t[counter]), str(t[counter])+"_cher", 400, 0., 200.)	
		RecEnergyHist = TH1F("RecEnergy_"+str(t[counter]),str(t[counter])+"_Energy", 400, 0., 200.)
		
		#Signalscinhist = TH1F("scintot_", str(counter+1)+"_scin", 3000, 0., 30000)
		EnergyHist = TH1F("Energy_"+str(t[counter]),str(t[counter])+"_Energy", 400, 0., 200.)
		LeakageHist = TH1F("Leak_"+str(t[counter]),str(t[counter])+"_Leak", 1000, 0., 100.)
		NeutrinoLeakageHist = TH1F("NeutrinoNeutrinoLeak_"+str(t[counter]),str(t[counter])+"_Leak", 1000, 0., 100.)
		TotalLeakageHist = TH1F("TotalLeak_"+str(t[counter]),str(t[counter])+"_Leak", 1000, 0., 100.)
		ChiHist = TH1F("Chi_"+str(t[counter]),str(t[counter])+"_Chi", 200, 0., 2.)
		scatterplot = TH2F("scatterplot_"+str(t[counter]), str(t[counter]), int(400), 0., 200., int(400), 0., 200.)
		EnergyContHist = TH1F("EnergyCont_"+str(t[counter]),str(t[counter])+"_EnergyCont", 400, 0., 200.)
		hesscatterplot = TH2F("H/E_S"+str(t[counter]), "H/E_S"+str(t[counter]), 200, 0., 1.1, 200, 0., 1.1)
		hecscatterplot = TH2F("H/E_C"+str(t[counter]), "H/E_C"+str(t[counter]), 200, 0., 1.1, 200, 0., 1.1)

		#loop over events
		entries = 50000
		for Event in range(entries):	

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
			
			#totalsignalscin = sum(BarrelR_VectorSignals)+sum(BarrelL_VectorSignals)
			#Signalscinhist.Fill(totalsignalscin)

			energytot = (sum(VectorR)+sum(VectorL))/1000
			EnergyHist.Fill(energytot)
			LeakageHist.Fill(Leak/1000.)
			NeutrinoLeakageHist.Fill(NeutrinoLeak/1000.)
			TotalLeakageHist.Fill(Leak/1000.+NeutrinoLeak/1000.)

			if float(t[counter])<12.:
				cutleak = 0.5
			if float(t[counter])>12. and float(t[counter])<50.:
				cutleak = 1.0
			if float(t[counter])>50.:
				cutleak = 3.0

			
			if (Leak/1000.+NeutrinoLeak/1000.)<cutleak:
				#apply calibrations
				Calib_BarrelL_VectorSignals = calibration.calibscin(BarrelL_VectorSignals)
				Calib_BarrelR_VectorSignals = calibration.calibscin(BarrelR_VectorSignals)
				Calib_BarrelL_VectorSignalsCher = calibration.calibcher(BarrelL_VectorSignalsCher)
				Calib_BarrelR_VectorSignalsCher = calibration.calibcher(BarrelR_VectorSignalsCher)
				#end of calibrations	
				energyscin = sum(Calib_BarrelR_VectorSignals)+sum(Calib_BarrelL_VectorSignals)
				energycher = sum(Calib_BarrelR_VectorSignalsCher)+sum(Calib_BarrelL_VectorSignalsCher)	
				e_c = float(t[counter])-(Leak/1000.+NeutrinoLeak/1000.)
				
				hesscatterplot.Fill(Energyem/1000./e_c, energyscin/e_c)
				hecscatterplot.Fill(Energyem/1000./e_c, energycher/e_c)

				EnergyContHist.Fill(e_c)
				ScinEnergyHist.Fill(energyscin)
				CherEnergyHist.Fill(energycher)
				scatterplot.Fill(energyscin/float(t[counter]), energycher/float(t[counter]))
				chi = chivalue
				newchi = (energyscin-e_c)/(energycher-e_c)
				ChiHist.Fill(newchi)			
				RecEnergyHist.Fill(1./0.99*(energyscin - chi*energycher)/(1.- chi))	
				
		print energies[counter], ScinEnergyHist.GetMean(), CherEnergyHist.GetMean(), RecEnergyHist.GetMean()
		displayfile.cd()
		gStyle.SetOptStat(111)
		#ScinEnergyHist.Fit("gaus")
		#CherEnergyHist.Fit("gaus")
		#RecEnergyHist.Scale(1/RecEnergyHist.Integral())
		#print RecEnergyHist.Integral()
		RecEnergyHist.Fit("gaus")
		RecEnergyHist.Write()
		ScinEnergyHist.Write()
		CherEnergyHist.Write()
		#Signalscinhist.Write()
		EnergyHist.Write()
		EnergyContHist.Write()
		e_cont = EnergyContHist.GetMean()
		e_cont_error = EnergyContHist.GetRMS()/float(entries)**0.5
		scatterplot.Write()
		#cut1 = TCutG("cut1",4)
		#cut1.SetVarX("x")
		#cut1.SetVarY("y")
		#cut1.SetPoint(0,0.,0.)
		#cut1.SetPoint(1,1.,0.)
		#cut1.SetPoint(2,1.,1.)
		#cut1.SetPoint(3,0.,1.)
		#profile = scatterplot.ProfileX("",1,400,"[cut1]")
		#profile.GetYaxis().SetRangeUser(0.,1.5)
		#profile.Write()
		func2 = TF1("func2", '[0]+x*(1.-[0])', 0., 1.)
		pp = hesscatterplot.ProfileX()
		pp.Fit(func2)
		pp.Write()
		hesscatterplot.Write()
		hecscatterplot.Write()
		ppc = hecscatterplot.ProfileX()
		ppc.Fit(func2)
		ppc.Write()
		LeakageHist.Write()
		NeutrinoLeakageHist.Write()
		TotalLeakageHist.Write()
		ChiHist.Write()
		#scin_sqrtenergies.append(1./(ScinEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		#cher_sqrtenergies.append(1./(CherEnergyHist.GetFunction("gaus").GetParameter(1)**0.5))
		Energy.append(RecEnergyHist.GetFunction("gaus").GetParameter(1))
		Energyerror.append(3.*RecEnergyHist.GetFunction("gaus").GetParameter(2)/float(entries)**0.5)
		MeanEnergyScin.append(ScinEnergyHist.GetMean())
		MeanEnergyCher.append(CherEnergyHist.GetMean())
		resolution.append(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))
		sigma_energy_e = ((RecEnergyHist.GetFunction("gaus").GetParameter(2)/(float(entries)**0.5))/RecEnergyHist.GetFunction("gaus").GetParameter(1) + RecEnergyHist.GetFunction("gaus").GetParError(2)/RecEnergyHist.GetFunction("gaus").GetParameter(2))*(RecEnergyHist.GetFunction("gaus").GetParameter(2)/RecEnergyHist.GetFunction("gaus").GetParameter(1))
		resolutionerror.append(sigma_energy_e)
		energyfractionscin.append(ScinEnergyHist.GetMean()/float(t[counter]))
		energyfractionscinerror.append(3*(ScinEnergyHist.GetRMS()/entries**0.5) / float(t[counter]))
		energyfractioncher.append(CherEnergyHist.GetMean()/ float(t[counter]))
		energyfractionchererror.append(3*(CherEnergyHist.GetRMS()/entries**0.5) / float(t[counter]))
		energyfraction.append(RecEnergyHist.GetFunction("gaus").GetParameter(1)/ float(t[counter]))
		energyfractionerror.append(3*(RecEnergyHist.GetFunction("gaus").GetParameter(2)/entries**0.5) / float(t[counter]))
		resolutionscin.append(ScinEnergyHist.GetRMS()/ScinEnergyHist.GetMean())
		sigma_energy_e = (ScinEnergyHist.GetMeanError()/ScinEnergyHist.GetMean() + ScinEnergyHist.GetRMSError()/ScinEnergyHist.GetRMS())*(ScinEnergyHist.GetRMS()/ScinEnergyHist.GetMean())
		resolutionscinerror.append(sigma_energy_e)
		resolutioncher.append(CherEnergyHist.GetRMS()/CherEnergyHist.GetMean())
		sigma_energy_e = (CherEnergyHist.GetMeanError()/CherEnergyHist.GetMean() + CherEnergyHist.GetRMSError()/CherEnergyHist.GetRMS())*(CherEnergyHist.GetRMS()/CherEnergyHist.GetMean())
		resolutionchererror.append(sigma_energy_e)
		chiarray.append(ChiHist.GetBinCenter(ChiHist.GetMaximumBin()))
		chierrorarray.append(3*ChiHist.GetRMS()/(entries**0.5))
		zeros.append(0.0)
		containment.append(e_cont/float(t[counter]))
		containmenterror.append(e_cont_error/float(t[counter]))

	containmentgraph = TGraphErrors(len(energies), energies, containment, zeros, containmenterror)
	containmentgraph.SetName("containment")
	containmentgraph.Write()
	ChiGraph = TGraphErrors(len(energies), energies, chiarray, zeros, chierrorarray)
	ChiGraph.SetName("Chi")
	ChiGraph.Write()
	LinearityGraph2 = TGraphErrors(len(energies), energies, Energy, zeros, Energyerror)
	LinearityGraph2.SetName("LinearityGraph2")
	LinearityGraph2.Write()
	LinearityGraph = TGraphErrors(len(energies), energies, energyfraction, zeros, energyfractionerror)
	LinearityGraph.SetName("LinearityGraph")
	LinearityGraph.Write()
	LinearityGraphScin = TGraphErrors(len(energies), energies, energyfractionscin, zeros, energyfractionscinerror)
	LinearityGraphCher = TGraphErrors(len(energies), energies, energyfractioncher, zeros, energyfractionchererror)
	LinearityGraphCher.SetName("LinearityGraphCher")
	LinearityGraphCher.Write()
	LinearityGraphScin.SetName("LinearityGraphScin")
	LinearityGraphScin.Write()

	ResolutionGraphScin = TGraphErrors(len(energies), sqrtenergies, resolutionscin, zeros, resolutionscinerror)
	func = TF1("func", "[0]/(x**0.5)+[1]", 10., 150.)
	ResolutionGraphCher = TGraphErrors(len(energies), sqrtenergies, resolutioncher, zeros, resolutionchererror)
	#ResolutionGraphScin.Fit("func", "R")
	#ResolutionGraphCher.Fit("func", "R")
	ResolutionGraphScin.SetName("ResolutionGraphScin")
	ResolutionGraphScin.Write()
	ResolutionGraphCher.SetName("ResolutionGraphCher")
	ResolutionGraphCher.Write()
	ResolutionGraph = TGraphErrors(len(energies), sqrtenergies, resolution, zeros, resolutionerror)
	#ResolutionGraph.Fit("func", "R")
	ResolutionGraph.SetName("ResolutionGraph")
	ResolutionGraph.Write()
	
	LinearityGraph.SetMinimum(0.976)#0.976
	LinearityGraph.SetMaximum(1.024)#1.024
	LinearityGraph.GetXaxis().SetLimits(0.0,155.0)
	LinearityGraph.SetTitle("")
	LinearityGraph.GetXaxis().SetLabelSize(.105)
	LinearityGraph.GetYaxis().SetLabelSize(0.105)
	LinearityGraph.GetXaxis().SetNdivisions(520)
	LinearityGraph.GetYaxis().SetNdivisions(504)
	LinearityGraph2.SetTitle("")
	LinearityGraph2.SetMinimum(0.0)
	LinearityGraph2.SetMaximum(155.0)
	LinearityGraph2.GetXaxis().SetLimits(0.0,155.0)
	LinearityGraph2.GetXaxis().SetLabelSize(0.)
	LinearityGraph2.GetXaxis().SetNdivisions(520)
	LinearityGraph2.GetYaxis().SetNdivisions(520)
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
	LinearityGraph2.Draw()
	pad2.cd()
	LinearityGraph.Draw("AP")
	ratioline = TF1("ratioline", str(np.mean(energyfraction)), 0., 160.)
	ratioline.SetLineColor(1)
	ratioline.SetLineWidth(1)
	ratioline.SetLineStyle(9)
	ratioline.Draw("same")
	ratioline.Write()
	c.Update()
	#c.SaveAs("MLratio.pdf")
	c.Write()


	'''	
	Resolutions = TMultiGraph()
	Resolutions.Add(ResolutionGraphScin)
	Resolutions.Add(ResolutionGraphCher)
	Resolutions.Add(ResolutionGraph)
	Resolutions.SetName("EMResolutions")
	Resolutions.Write()

	Linearities = TMultiGraph()
	Linearities.Add(LinearityGraph)
	Linearities.Add(LinearityGraphScin)
	Linearities.Add(LinearityGraphCher)
	Linearities.SetName("Linearities")
	Linearities.Write()
	'''	
	return RecEnergyHist.GetFunction("gaus").GetParameter(2), RecEnergyHist.GetFunction("gaus").GetParameter(1), RecEnergyHist.GetFunction("gaus").GetChisquare()
#names = ["FTFPBERT"]#, "FTFPBERT", "QGSPBERT", "QBBC"]
#names = ["FTFPBERTTRV","FTFPBERT", "QGSPBERT", "QBBC"]
name = raw_input("physics list: ")
names = []
names.append(name)
for name in names:
	recenergy(name, 0.41)
'''
#for results varying chi
for name in names:
	chis = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7]
	sigmas = array('d')
	energies = array('d')
	sigma_energies = array('d')
	chisquare_array = array('d')
	for chi in chis:
		print chi
		sigma, energy, chisquare = recenergy(name, chi)
		sigmas.append(sigma)
		energies.append(energy)
		sigma_energies.append(sigma/energy)
		chisquare_array.append(chisquare)
	outfile = TFile("newout.root","RECREATE")
	TGraph1 = TGraph(int(len(chis)), array('d',chis),sigmas)
	TGraph1.SetName("sigmas")
	TGraph1.SetTitle("sigmas")
	TGraph1.Write()
	TGraph2 = TGraph(int(len(chis)), array('d',chis),energies)
	TGraph2.SetName("energies")
	TGraph2.SetTitle("energies")
	TGraph2.Write()
	TGraph3 = TGraph(int(len(chis)), array('d',chis),sigma_energies)
	TGraph3.SetTitle("sigma_energy")
	TGraph3.Write()
	TGraph4 = TGraph(int(len(chis)), array('d',chis), chisquare_array)
	TGraph4.SetTitle("chisquare")
	TGraph4.Write()
'''