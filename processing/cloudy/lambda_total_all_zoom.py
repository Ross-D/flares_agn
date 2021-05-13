import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import matplotlib as mpl
import matplotlib.cm as cm

output_dir = 'output/Log'

AGN_T = np.logspace(4, np.log10(1.5e6), 100)

cmap = mpl.cm.viridis
norm = mpl.colors.LogNorm(vmin = 1e4, vmax = 1.5e6)

fig = plt.figure(figsize = (14,8))
for i in range(100):
    i+=1
    cont = np.loadtxt(f'{output_dir}/{i}.cont', delimiter='\t', dtype=str)#, usecols = (0,1,2,3,4,8))
    cont_nu, cont_incident, cont_transmitted, cont_nebular, cont_total, cont_linecont = cont[:,0], cont[:,1], cont[:,2], cont[:,3], cont[:,4], cont[:,8], 

    nu, incident, transmitted, total = np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu))
    for ii in range(len(nu)):
        nu[ii] = float(cont_nu[ii])
        incident[ii] = float(cont_incident[ii])
        transmitted[ii] = float(cont_transmitted[ii])
        total[ii] = float(cont_total[ii])
    
    plt.plot(nu, total, c = cmap(norm(AGN_T[i-1])), label = str(i))
plt.yscale('log')
plt.xscale('log')
#plt.ylim(1e37, 1e45)
plt.xlim(1e1, 1e6)
plt.title("Total Spectrum")
plt.ylabel(r"$\lambda F_\lambda$ / ($10^{42}$ erg s$^-1$)")
plt.xlabel(r"$\lambda$ / $\AA$")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), label = "Temperature /K")#, ax=ax)
#plt.legend()
fig.savefig("Plots/Log/lambda_total_all_zoom.png")
fig.clf()

for i in range(100):
    i+=1
    cont = np.loadtxt(f'{output_dir}/{i}.cont', delimiter='\t', dtype=str)#, usecols = (0,1,2,3,4,8))
    cont_nu, cont_incident, cont_transmitted, cont_nebular, cont_total, cont_linecont = cont[:,0], cont[:,1], cont[:,2], cont[:,3], cont[:,4], cont[:,8], 

    nu, incident, transmitted, total = np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu))
    for ii in range(len(nu)):
        nu[ii] = float(cont_nu[ii])
        incident[ii] = float(cont_incident[ii])
        transmitted[ii] = float(cont_transmitted[ii])
        total[ii] = float(cont_total[ii])
    
    plt.plot(nu, total, c = cmap(norm(AGN_T[i-1])), label = str(i))
plt.yscale('log')
plt.xscale('log')
#plt.ylim(1e37, 1e45)
plt.xlim(1e1, 9e2)
plt.title("Total Spectrum")
plt.ylabel(r"$\lambda F_\lambda$ / ($10^{42}$ erg s$^-1$)")
plt.xlabel(r"$\lambda$ / $\AA$")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), label = "Temperature /K")#, ax=ax)
#plt.legend()
fig.savefig("Plots/Log/lambda_total_all_zoom_low.png")
fig.clf()

for i in range(100):
    i+=1
    cont = np.loadtxt(f'{output_dir}/{i}.cont', delimiter='\t', dtype=str)#, usecols = (0,1,2,3,4,8))
    cont_nu, cont_incident, cont_transmitted, cont_nebular, cont_total, cont_linecont = cont[:,0], cont[:,1], cont[:,2], cont[:,3], cont[:,4], cont[:,8], 

    nu, incident, transmitted, total = np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu))
    for ii in range(len(nu)):
        nu[ii] = float(cont_nu[ii])
        incident[ii] = float(cont_incident[ii])
        transmitted[ii] = float(cont_transmitted[ii])
        total[ii] = float(cont_total[ii])
    
    plt.plot(nu, total, c = cmap(norm(AGN_T[i-1])), label = str(i))
plt.yscale('log')
plt.xscale('log')
plt.ylim(1e33, 0.5e46)
plt.xlim(8e2, 5e6)
plt.title("Total Spectrum")
plt.ylabel(r"$\lambda F_\lambda$ / ($10^{42}$ erg s$^-1$)")
plt.xlabel(r"$\lambda$ / $\AA$")
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), label = "Temperature /K")#, ax=ax)
#plt.legend()
fig.savefig("Plots/Log/lambda_total_all_zoom_high.png")
fig.clf()
