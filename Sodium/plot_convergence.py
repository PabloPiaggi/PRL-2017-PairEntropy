import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

###################################################################
# Plot options
###################################################################
font = {'family' : 'serif',
        'serif'   : 'palatino',
        'style'   : 'normal',
        'variant'   : 'normal',
        'stretch'   : 'normal',
        'weight'   : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = (8.3, 5.8)

###################################################################
# Colors
###################################################################
# Palette 1
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
MarkerSize=20
palette=tableau20

###################################################################
# Plot data
###################################################################

temp=[300,325,350,375,400]
fileList=['300K-MetaD/Convergence/results.txt','325K-MetaD/Convergence/results.txt','350K-MetaD/Convergence/results.txt','375K-MetaD/Convergence/results.txt','400K-MetaD/Convergence/results.txt']

counter=0
for file in fileList:
	deltaG = np.loadtxt(file, unpack=True)
	mean = np.mean(deltaG[50:])
	std = np.std(deltaG[50:])
	print(mean,std)
	time = np.linspace(0,deltaG.shape[0],deltaG.shape[0])+1
	plt.plot([np.amin(time),np.amax(time)],[mean,mean],'--',color="black",linewidth=1.0,alpha=0.6)
	plt.fill_between([np.amin(time),np.amax(time)],mean-std,mean+std,color=palette[counter*2],linewidth=1.0,alpha=0.5)
	plt.plot(time,deltaG,'-',color=palette[counter*2],linewidth=3.0,alpha=0.8)
	plt.text(400, mean-std-4, str(temp[counter])+ " K",va='center',ha='center',color=palette[counter*2])
	counter+=1

###################################################################
# Plot options
###################################################################

plt.xlim([0,500])
plt.ylim([-80,140])
plt.xlabel("Time (ns)")
plt.ylabel(r'$\Delta$$G_{S \rightarrow L}$ (kJ/mol)')
plt.tick_params(axis='x', pad=10)
plt.tick_params(axis='y', pad=10)
plt.xticks(np.linspace(0,500,6),np.linspace(0,500,6).astype(int))
plt.show()
