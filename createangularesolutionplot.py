from ROOT import TF1, TGraph, TFile, gStyle
from array import array

energies_arr = array('d')
angrestheta_arr = array('d')
angresphi_arr = array('d')

energies = [1.,10., 20., 30., 40., 50., 70., 80., 90., 110., 130., 150.]
angrestheta = [3.582, 2.884, 2.175, 1.739, 1.581, 1.442, 1.225, 1.1112, 1.064, 0.8892, 0.7539, 0.7463]
angresphi = [2.811, 1.275, 0.8831, 0.7306, 0.6509, 0.5869, 0.4639, 0.4436, 0.3966, 0.3573, 0.3128, 0.2962]

outputfile = TFile("angres.root", "RECREATE")
for counter, i in enumerate(angresphi):
	angresphi_arr.append(i)
for counter, t in enumerate(angrestheta):
	angrestheta_arr.append(t)
for s in energies:
	energies_arr.append(s)
print angrestheta_arr
print angresphi_arr
ThetaGraph = TGraph(12, energies_arr, angrestheta_arr)
PhiGraph = TGraph(12, energies_arr, angresphi_arr)

resolutionfittheta = TF1("resolutionfit", '[0]/(x**0.5)+[1]/x+[2]', 0.5, 160.)
resolutionfittheta.SetParLimits(2, 0.5, 0.8)
resolutionfitphi = TF1("resolutionfit", '[0]/(x**0.5)+[1]/x+[2]', 0.5, 160.)
resolutionfitphi.SetParLimits(2, 0.5, 0.8)

ThetaGraph.Fit(resolutionfittheta)
PhiGraph.Fit(resolutionfitphi)

gStyle.SetOptStat(1)
ThetaGraph.Write()
PhiGraph.Write()



