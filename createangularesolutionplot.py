from ROOT import TF1, TGraph, TFile, gStyle, TGraphErrors
from array import array

#energies = [1.,10., 20., 30., 40., 50., 70., 80., 90., 110., 130., 150.]
#angrestheta = [3.582, 2.884, 2.175, 1.739, 1.581, 1.442, 1.225, 1.1112, 1.064, 0.8892, 0.7539, 0.7463]
#angresphi = [2.811, 1.275, 0.8831, 0.7306, 0.6509, 0.5869, 0.4639, 0.4436, 0.3966, 0.3573, 0.3128, 0.2962]

def angresplot(energies, angrestheta, angresphi, angresthetaerror, angresphierror):

	energies_arr = array('d', energies)
	angrestheta_arr = array('d', angrestheta)
	angresphi_arr = array('d', angresphi)
	angresthetaerror_arr = array('d', angresthetaerror)
	angresphierror_arr = array('d', angresphierror)
	zeros_arr = array('d',[0.0]*len(energies))

	ThetaGraph = TGraphErrors(len(energies), energies_arr, angrestheta_arr, zeros_arr, angresthetaerror_arr)
	PhiGraph = TGraphErrors(len(energies), energies_arr, angresphi_arr, zeros_arr, angresphierror_arr)	

	resolutionfittheta = TF1("resolutionfit", '[0]/(x**0.5)+[1]', 30., 150.)
	#resolutionfittheta.SetParLimits(2, 0.5, 0.8)
	resolutionfitphi = TF1("resolutionfit", '[0]/(x**0.5)+[1]', 30., 150.)
	#resolutionfitphi.SetParLimits(2, 0.5, 0.8)	

	ThetaGraph.Fit(resolutionfittheta)
	PhiGraph.Fit(resolutionfitphi)	
	return ThetaGraph, PhiGraph



