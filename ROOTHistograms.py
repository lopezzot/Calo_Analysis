from ROOT import *
import map
import newmap

def create_eventdisplay(PrimaryParticleName, VectorSignalsR, VectorSignalsL, filename):
	"""Function to perform ROOT event display form calo"""
	NbOfBarrel=40
	NbOfEndcap=35
	NZrot=36
	TotTower=NbOfBarrel+NbOfEndcap

	deltatheta = 45./(NbOfBarrel)

	#Set ROOT histograms (x=theta, y=phi)
	TH2Signals = TH2F("ScinSignals"+filename, PrimaryParticleName, NbOfBarrel*2, -1*(deltatheta*NbOfBarrel+deltatheta/2), deltatheta*NbOfBarrel+deltatheta/2, NZrot, 0., 360.)
	
	#Fill histograms in for loop
	for towerindex in range(TotTower*NZrot):
		theta, phi = newmap.maptower(towerindex, "right")
		TH2Signals.Fill(theta,phi, VectorSignalsR[towerindex])
		theta, phi = newmap.maptower(towerindex, "left")
		TH2Signals.Fill(theta, phi, VectorSignalsL[towerindex])

	#Draw + DrawOptions histograms	

	Style = gStyle
	Style.SetPadRightMargin(0.16)
	Style.SetPalette(1) #Root palette style
	Style.SetOptStat(0) #Do not show statistics
	TH2Signals.SetLineWidth(0) #TH2Signals #No line width
	TH2Signals.SetLineColor(2)
	#TH2Signals.SetFillColorAlpha(2, 0.)
	XAxis = TH2Signals.GetXaxis()
	XAxis.SetTitle("Theta (deg)")
	XAxis.CenterTitle()
	#XAxis.SetTitleOffset(1.8)
	YAxis = TH2Signals.GetYaxis()
	YAxis.SetTitle("Phi (deg)")
	YAxis.CenterTitle()
	#YAxis.SetTitleOffset(1.8)
	ZAxis = TH2Signals.GetZaxis()
	#ZAxis.SetTitle("Energy (MeV)")
	#ZAxis.SetTitleOffset(1.4)
	TH2Signals.Draw("COLZ 0 FB")
	TH2Signals.Write()
	gPad.SaveAs("ImageScintillation"+filename+".pdf")
	