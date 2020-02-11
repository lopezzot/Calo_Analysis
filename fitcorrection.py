from ROOT import *
from array import array

#sigma = array('d', (1.517, 1.914, 2.237, 2.544))
#rms = array('d', (1.9, 2.45, 2.92, 3.15, 5.63))
#e = array('d', (30, 50, 70, 90, 250))
#y = array('d', [s/e[counter] for counter, s in enumerate(rms)])
y = array('d', (0.3987,0.6025,0.8405,1.09))
y = array('d', (0.7836, 1.211, 1.65, 2.106)) #fit for chi = 0.29
y = array('d', (0.54, 0.67, 0.75, 0.84)) #fit for X0 correction
y = array('d', (0.7153, 1.071, 1.45, 1.754, 2.634, 3.85)) #fit for X0 correction
y = array('d', (0.73, 1.13, 1.56, 1.93, 3.11, 5.3)) #fit for X0 correction
y = array('d', (1.237, 1.734, 2.163, 2.586, 3.615, 5.039)) #fit for X0 correction
y = array('d', (1.26, 1.81, 2.31, 2.82, 4.25, 6.7)) #fit for X0 correction

x = array('d', (15.,25.,35.,45., 75.,125.))

fitfile = TFile("fitcorrection.root", "RECREATE")
graph = TGraph(6, x, y)

func = TF1("func", "[0]*x+[1]", 15., 125.)
graph.Fit("func", "R")
graph.Write()
	