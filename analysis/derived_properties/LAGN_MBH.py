def L_AGN(Mdot):
    # This function will calculate the bolometric luminosity of an
    # AGN from an input of accretion rate input of log10(Mdot/(M_sol/yr))
    # the output will be log10(Lagn/(erg/s))
    import numpy as np
    import astropy.constants as const
    import astropy.units as u
    Msun = const.M_sun
    Lsun = const.L_sun
    c = const.c
    Mdot_SI = ((10**Mdot)*Msun*(u.yr**-1))
    n = 0.1
    Lagn_SI = n*Mdot_SI*(c**2)
    Lagn_CGS = Lagn_SI.to(u.erg/u.s)
    Lagn = np.log10(Lagn_CGS/(u.erg/u.s))
    return Lagn

def L_edd(M_BH):
    # This function will calculate the Eddington luminosity of an
    # AGN from an input of accretion rate input of log10(MBH/(M_sol))
    # the output will be log10(Ledd/(erg/s))
    import numpy as np
    import astropy.constants as c
    import astropy.units as u
    MBH_= ((10**M_BH)*u.M_sun)
    MBH_SI = MBH_.to(u.kg)
    numerator = 4*np.pi*c.G*c.m_p*c.c*MBH_SI
    Ledd = numerator/c.sigma_T
    L_edd = Ledd.to(u.erg/u.s)
    return (L_edd/(u.erg/u.s))
    
def Lagn_MBH_grid():
    '''No inputs, creates a figure of 6 plots of Lagn vs MBH at each redshift'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units

    flares_dir = '../../../../../../../research/astro/flare' #location of flares directory (apollo2 location specified)
    
    zed_tags = (-1,-2,-3,-4,-5,-6) #labels used to specify redshifts in flares
    
    #importing flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #specifying colormap
    cmap = mpl.cm.plasma
    
    #loading FLARES data
    Mstar = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy') # Black accretion rate
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy

    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    temp_max = []
    temp_min = []
    
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag] 
        z = fl.zeds[zed_tag] #converts zed tag to redshift

        #getting weights from flares
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        #creating data arrarys
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(Mstar[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        
        MBH += 10. # convert MBH units to M_sol
        h = 0.6777 # Hubble parameter
        
        #converting MBHacc units to M_sol/yr
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr

        x += 10 # convert Mstar units to M_sol
        
        #cutting out galaxies below Mstar = 10**8
        y = MBH
        y_dot = MBHacc
        
        xcut, ycut, y_dotcut = np.array([]), np.array([]), np.array([])
        for i in range(len(y)):
            if x[i] > 8:
                xcut = np.append(xcut, x[i])
                ycut = np.append(ycut, y[i])
                y_dotcut = np.append(y_dotcut, y_dot[i])
            else:
                continue
        
        x = xcut
        y = ycut
        y_dot = y_dotcut
        
        #computing luminosities with function
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
         
        #Creating plot for this redshift
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        MBHrange = np.linspace(5.5, 10, 100)
        ax.scatter(lums, y, alpha=0.25)
        ax.plot(np.log10(L_edd(MBHrange)), MBHrange, c = 'k', linestyle = 'dashed')
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/erg s^{-1}))$')
        ax.set_ylabel(r'SMBH Mass ($\rm log_{10}(M_{BH}/M_{\odot}))$')
        ax.set_xlim(35,48)
        ax.set_ylim(5.5,10)
        ax.grid()
        plot_no+=1

    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/LAGN_MBH/LAGN_MBH_grid_erg.png')
    plt.clf()
    fig.clf()
    
def Lagn_MBH_individual():
    '''No inputs, creates an individual plot of Lagn vs MBH at each redshift'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units

    flares_dir = '../../../../../../../research/astro/flare'    # location of flares directory (apollo2 location specified)
    
    zed_tags = (-1,-2,-3,-4,-5,-6) # labels used to specify redshifts in flares
    
    # importing datasets from flares
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #specifying colormap
    cmap = mpl.cm.plasma
    
    #loading datasets from flares
    Mstar = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy') # Black accretion rate
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy

    X = Mstar
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    temp_max = []
    temp_min = []
    
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag] # convertin zed tag to redshift number
        
        #loading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        #creating data arrays
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        MBH += 10. # converting MBH to M_sol

        # converting MBHacc to M_sol/yr
        h = 0.6777 # Hubble parameter
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr

        x += 10 # converting Mstar to M_sol
        
        #cutting out galaxies below mstar = 10**8
        y = MBH
        y_dot = MBHacc
        
        xcut, ycut, y_dotcut = np.array([]), np.array([]), np.array([])
        for i in range(len(y)):
            if x[i] > 8:
                xcut = np.append(xcut, x[i])
                ycut = np.append(ycut, y[i])
                y_dotcut = np.append(y_dotcut, y_dot[i])
            else:
                continue
        
        x = xcut
        y = ycut
        y_dot = y_dotcut
        
        #calculating luminosities with function
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
        
        #plotting Lagn against MBH
        MBHrange = np.linspace(5, 10, 100)
        plt.scatter(lums, y, alpha=0.25)
        plt.plot(np.log10(L_edd(MBHrange)), MBHrange, c = 'k', linestyle = 'dashed')
        plt.title('z = '+str(z))
        plt.xlabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/erg s^{-1}))$')
        plt.ylabel(r'SMBH Mass ($\rm log_{10}(M_{BH}/M_{\odot}))$')
        plt.xlim(35)
        plt.ylim(5.16,10)
        plt.grid()
        plot_no+=1

        plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
        fig.subplots_adjust(right=0.85)
        fig.savefig(f'figures/LAGN_MBH/LAGN_MBH_'+str(int(z))+'_erg.png')
        plt.clf()
        fig.clf()
