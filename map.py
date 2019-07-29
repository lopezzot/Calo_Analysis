import math

def maptowerBR(index):
	"""Function to return tower number given index in BarrelRvectorsignals"""
	sliceindex = index//40
	towerindex = index-(sliceindex*40)

	deltatheta_barrel = [0.02222,0.02220,0.02217,0.02214,0.02209,0.02203,0.02196,0.02188,0.02179,0.02169,0.02158,0.02146,0.02133,0.02119,0.02105,0.02089,0.02073,0.02056,0.02039,0.02020,0.02002,0.01982,0.01962,0.01941,0.01920,0.01898,0.01876,0.01854,0.01831,0.01808,0.01785,0.01761,0.01738,0.01714,0.01689,0.01665,0.01641,0.01616,0.01592,0.01567]
	deltatheta_barrel = [math.degrees(x) for x in deltatheta_barrel]

	#get theta
	theta = sum(deltatheta_barrel[0:towerindex])

	#get phi
	phi_unit = math.degrees(2.*math.pi/256.) #256 is number of Zrot 
	phi = sliceindex*phi_unit

	return theta, phi