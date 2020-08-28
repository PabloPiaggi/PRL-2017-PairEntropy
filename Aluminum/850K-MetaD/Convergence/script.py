import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import cm
#from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate

###################################################################
# Set constants
###################################################################
plotMaxValue=500.
watershed=-309.0

###################################################################
# Calculation
###################################################################

fileBase='fes_'
temperature=850
numParticles=256
deltaG=[]
for fileNumber in range(1,1001,1):
	x, y, z, dummy1, dummy2 = np.loadtxt(fileBase + str(fileNumber) + ".dat", unpack=True)
	N = int(len(z)**.5)
	x=x[0:N]
	energyDivision=np.argmin(np.abs(x-watershed))
	z = z.reshape(N, N)
	z = np.minimum(z, plotMaxValue)
	beta=1./(0.0083144621*temperature)
	z=np.exp(-beta*z)
	integral=np.trapz(z,axis=0) # integrate in entropy first
	#plt.plot(-(1./beta)*np.log(integral))
	# Solid Basin
	integral2=np.trapz(integral[:energyDivision])
	freeEnergySolid=-(1./beta)*np.log(integral2)
	# Liquid Basin
	integral3=np.trapz(integral[energyDivision:])
	freeEnergyLiquid=-(1./beta)*np.log(integral3)
	#print(freeEnergySolid,freeEnergyLiquid,freeEnergyLiquid-freeEnergySolid)
	deltaG += freeEnergyLiquid-freeEnergySolid,
#plt.plot([energyDivision,energyDivision],[0,500],'--')
np.savetxt('results.txt',deltaG)
#plt.show()

