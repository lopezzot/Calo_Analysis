from ROOT import TF1, TGraph, TFile
from array import array

energies_arr = array('d')
angrestheta_arr = array('d')
angresphi_arr = array('d')

energies = [1.,10., 40., 80., 100.]
angrestheta = [1.523, 0.709, 0.4641, 0.1988, 0.1859]
angresphi = [1.113, 0.3696, 0.2016, 0.1141, 0.1146]
outputfile = TFile("angres.root", "RECREATE")
for counter, i in enumerate(angresphi):
	angresphi_arr.append(i)
for counter, t in enumerate(angrestheta):
	angrestheta_arr.append(t)
for s in energies:
	energies_arr.append(s)
print angrestheta_arr
print angresphi_arr
ThetaGraph = TGraph(5, energies_arr, angrestheta_arr)
PhiGraph = TGraph(5, energies_arr, angresphi_arr)

resolutionfittheta = TF1("resolutionfit", '[0]/(x**0.5)+[1]', 1., 100.)
resolutionfitphi = TF1("resolutionfit", '[0]/(x**0.5)+[1]', 1., 100.)

ThetaGraph.Fit(resolutionfittheta)
PhiGraph.Fit(resolutionfitphi)

ThetaGraph.Write()
PhiGraph.Write()



