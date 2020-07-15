from ROOT import *
import numpy as np
import math

def cart2sph(x,y,z):
	if y >= 0.:
		phi = np.arccos(x/np.sqrt(x**2+y**2))
	else:
		phi = 2*math.pi-np.arccos(x/np.sqrt(x**2+y**2))
	#phi = np.arctan2(y,x)
	theta = np.arctan2(z,np.sqrt(x**2 + y**2))
	r = np.sqrt(x**2 + y**2 + z**2)
	return theta, phi, r

Theta = 30
Phi = 90

Theta_rad = Theta*math.pi/180.
Phi_rad = Phi*math.pi/180.

X = math.sin(math.pi/2-Theta_rad)*math.cos(Phi_rad)
Y = math.sin(math.pi/2-Theta_rad)*math.sin(Phi_rad)
Z = math.cos(math.pi/2-Theta_rad)

print Theta, Phi
print Theta_rad, Phi_rad
print X,Y,Z
print cart2sph(X,Y,Z)
print cart2sph(X,Y,Z)[0]*180./math.pi, cart2sph(X,Y,Z)[1]*180./math.pi

#Best way of doing angles using root
Theta = 110
Phi = 350

Theta_rad = Theta*math.pi/180.
Phi_rad = Phi*math.pi/180.
a = TVector3()
a.SetMagThetaPhi(1.,Theta_rad,Phi_rad)
print a.X(),a.Y(),a.Z()
#end
