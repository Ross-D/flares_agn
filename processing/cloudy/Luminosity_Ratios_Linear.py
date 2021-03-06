import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.integrate import simps

import astropy.units as u
import astropy.constants as c
lower_lim = (c.c/(((1 * u.Ry).to(u.J))/c.h)).to(u.AA)
upper_lim = (c.c/((((7.354e6) * u.Ry).to(u.J))/c.h)).to(u.AA)

output_dir = 'output/Linear'

cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin = 1, vmax = 100)

AGN_T_ = np.linspace(1e4, 1.5e6, 100)
Integrals = {}
L_FUV = []
L_Tot = []
L_Inc = []
L_Trans = []
AGN_T = []
for i in range(100):
    
    temp = AGN_T_[i]
    temp = round(int(temp), 3 - len(str(int(temp))))
    i+=1
    cont = np.loadtxt(f'{output_dir}/{i}.cont', delimiter='\t', dtype=str)#, usecols = (0,1,2,3,4,8))
    cont_nu, cont_incident, cont_transmitted, cont_nebular, cont_total, cont_linecont = cont[:,0], cont[:,1], cont[:,2], cont[:,3], cont[:,4], cont[:,8], 

    nu, incident, transmitted, total = np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu)), np.zeros(len(cont_nu))
    for ii in range(len(nu)):
        nu[ii] = float(cont_nu[ii])
        incident[ii] = float(cont_incident[ii])
        transmitted[ii] = float(cont_transmitted[ii])
        total[ii] = float(cont_total[ii])
        
        
    l_tot, l_inc, l_trans, l_fuv, nu_fuv = [], [], [], [], []
    for i in range(len(nu)):
        l_tot.append(total[i]/nu[i])
        l_inc.append(incident[i]/nu[i])
        l_trans.append(transmitted[i]/nu[i])
        if nu[i] < upper_lim/u.AA:
            continue
        elif nu[i] > lower_lim/u.AA:
            continue
        else:
            nu_fuv.append(nu[i])
            l_fuv.append(total[i]/nu[i])
    '''print(('-'*10)+str(temp)+('-'*10))
    print(len(nu))
    print(len(nu_fuv))'''
    l_tot_int = simps(np.flip(l_tot), np.flip(nu))
    l_inc_int = simps(np.flip(l_inc), np.flip(nu))
    l_trans_int = simps(np.flip(l_trans), np.flip(nu))
    l_fuv_int = simps(np.flip(l_fuv), np.flip(nu_fuv))
    
    L_FUV.append(l_fuv_int)
    L_Tot.append(l_tot_int)
    L_Inc.append(l_inc_int)
    L_Trans.append(l_trans_int)
    AGN_T.append(temp)
        
Integrals.update({"L_FUV":L_FUV})
Integrals.update({"L_Tot":L_Tot})
Integrals.update({"L_Inc":L_Inc})
Integrals.update({"L_Trans":L_Trans})
Integrals.update({"AGN_T":AGN_T})

y1, y2 = [], []
for ii in range(len(Integrals["L_FUV"])):
    y1.append(Integrals["L_FUV"][ii]/Integrals["L_Inc"][ii])
    y2.append(Integrals["L_FUV"][ii]/Integrals["L_Tot"][ii])

x = []
for j in range(len(Integrals["AGN_T"])):
    x.append(Integrals["AGN_T"][j]/1e5)
    
plt.figure(figsize = (14,8))
plt.plot(x, y1)
plt.ylabel(r"$L_{FUV}/L_{Inc}$")
plt.xlabel(r"T$_{AGN}$/($10^5$ K)")
plt.title(r"Ratio of Observed far UV Luminosity to Incident Bolometric Luminosity")
plt.xlim(0,15)
plt.ylim(0.04,0.12)
plt.grid()
plt.savefig("Plots/Ratios/L_INC_L_TOT_Linear.png")

plt.figure(figsize = (14,8))
plt.plot(x, y2)
plt.ylabel(r"$L_{FUV}/L_{Bol}$")
plt.xlabel(r"T$_{AGN}$/($10^5$ K)")
plt.title(r"Ratio of Observed far UV and Bolometric Luminosities")
plt.xlim(0,15)
plt.ylim(0.13, 0.145)
plt.grid()
plt.savefig("Plots/Ratios/L_FUV_L_TOT_Linear.png")

def f(x_val, m_val, c_val):
    return m_val*(x_val)+c_val
from scipy.optimize import curve_fit as fit

x_fit = []
y_fit = []
for val in range(len(x)):
    if x[val] > 2.2 :
        x_fit.append(x[val])
        y_fit.append(y2[val])
    else:
        continue

popt, pcov = fit(f, x_fit, y_fit)



m = popt[0]
c = popt[1]

Fitted = []
Rounded_Fit = []
for i in range(len(x_fit)):
    Fitted.append(f(x_fit[i], m, c))
    Rounded_Fit.append(f(x_fit[i], -0.0006977, 0.1448))

def curve(x, a, b, c, d, e, f, g):
    return g*((a*(x**b)) + (c/(d+(x**e))) + f)
popt_T, pcov_T = fit(curve, x, y2)

y_T = []
for j in range(len(x)):
    y_T.append(curve(x[j], popt_T[0], popt_T[1], popt_T[2], popt_T[3], popt_T[4], popt_T[5], popt_T[6]))\

plt.figure(figsize = (14,8))
plt.plot(x, y2, "--", label = "Data")
plt.plot(x_fit, Fitted, label = r"Linear fit")
#plt.plot(x_fit, Rounded_Fit, "o", label = r"Rounded Linear fit")
plt.plot(x, y_T, label = "Full fit")
plt.ylabel(r"$L_{FUV}/L_{Bol}$")
plt.xlabel(r"T$_{AGN}$/($10^5$ K)")
plt.title(r"Ratio of Observed far UV and Bolometric Luminosities")
plt.ylim(0.13, 0.145)
plt.xlim(0, 15)
plt.legend()
plt.grid()
plt.savefig("Plots/Ratios/L_FUV_L_TOT_Linear_withfit.png")

print(popt)
