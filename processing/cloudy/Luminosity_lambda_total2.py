import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import matplotlib as mpl
import matplotlib.cm as cm

output_dir = 'output/Luminosity'

AGN_L = np.arange(40, 49)

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin = 40, vmax = 49)

#runs = 10*np.arange(1,10)
#runs = np.append(1, 10*np.arange(1,11,1))
fig = plt.figure(figsize = (14,8))
#for i in runs:
for i in range(9):
    i+=1
    cont = np.loadtxt(f'{output_dir}/{i}.cont', delimiter='\t', dtype=str)#, usecols = (0,1,2,3,4,8))
    cont_nu, cont_incident, cont_transmitted, cont_nebular, cont_total, cont_linecont = cont[:,0], cont[:,1], cont[:,2], cont[:,3], cont[:,4], cont[:,8], 

    nu, incident, transmitted, total = np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu))
    for ii in range(len(nu)):
        nu[ii] = float(cont_nu[ii])
        incident[ii] = float(cont_incident[ii])
        transmitted[ii] = float(cont_transmitted[ii])
        total[ii] = float(cont_total[ii])
    
    plt.plot(nu, total, c = cmap(norm(AGN_L[i-1])), label = "1e"+str(AGN_L[i-1])+r" $erg s^{-1}$")
plt.yscale('log')
plt.xscale('log')
#plt.ylim(1e33, 1e49)
plt.xlim(1e-2, 1e7)
plt.title("Total Spectrum")
plt.ylabel(r"$\lambda F_\lambda$ / ($10^{42}$ erg s$^-1$)")
plt.xlabel(r"$\lambda$ / $\AA$")
plt.legend()
fig.savefig("Plots/Luminosity/lambda_total2.png")
fig.clf()
