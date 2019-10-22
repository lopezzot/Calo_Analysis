import math
import numpy

def maptower(index, side):
	"""Function to return tower angles (theta and phi) given index"""
	NbOfBarrel=40
	NbOfEndcap=35
	NZrot=36
	TotTower=NbOfBarrel+NbOfEndcap
	index = index-1
	sliceindex = index//TotTower
	towerindex = index-(sliceindex*TotTower)
	deltatheta = 45./(NbOfBarrel)
	
	#get theta
	theta = towerindex*deltatheta+deltatheta/2.

	#get phi
	phi_unit = 360./NZrot
	
	phi = (sliceindex)*phi_unit

	'''
	if phi<180.:
		phi = 180+phi
	elif phi>180.:
		phi = phi-180
    '''
	if side == "right":
		#print theta+90., phi, -numpy.log(numpy.tan(((90.-theta)*math.pi/180./2.)))
		return theta+90., phi, -numpy.log(numpy.tan(((90.-theta)*math.pi/180./2.)))
	if side == "left":
		#print 90.-theta, phi, numpy.log(numpy.tan(((90.-theta)*math.pi/180./2.)))
		return 90.-theta, phi, numpy.log(numpy.tan(((90.-theta)*math.pi/180./2.)))

