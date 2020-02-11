from ROOT import *
from array import array

sigma = array('d', (1.5, 1.88, 2.283, 2.47))
#rms = array('d', (1.9, 2.45, 2.92, 3.15, 5.63))
e = array('d', (15., 25., 35., 45.))
y = array('d', [s/e[counter] for counter, s in enumerate(sigma)])


fitfile = TFile("fitfile.root", "RECREATE")
graph = TGraph(4, e, y)

func = TF1("func", "[0]/(x**0.5)", 15., 45.)
graph.Fit("func", "R")
graph.Write()
	