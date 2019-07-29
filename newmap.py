import math

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
	
	if phi<180.:
		phi = 180+phi
	elif phi>180.:
		phi = phi-180
    
	if side == "right":
		return theta, phi
	if side == "left":
		return -1*theta, phi