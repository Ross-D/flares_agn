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
    Lagn_SI = Mdot_SI*(c**2)
    Lagn = np.log10(Lagn_SI/Lsun)
    return Lagn

def Lagn_Mstar_individual(wanted_zed_tag):
    '''takes single input, must be string:
    'all' gives Lagn vs M* for all possible redshifts
    '-1' gives redshift 5
    '-2' gives redshift 6
    '-3' gives redshift 7
    '-4' gives redshift 8
    '-5' gives redshift 9
    '-6' gives redshift 10
        '''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units
    
    #sorting which z tag is wanted
    if wanted_zed_tag == 'all':
        zed_tags = (-1,-2,-3,-4,-5,-6)
    else:
        zed_tags = [int(wanted_zed_tag)]
        
    flares_dir = '../../../../../../../research/astro/flare' #location of flares
    
    # importing flars data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #specifying colormap
    cmap = mpl.cm.plasma
    
    #loading dataset
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
        z = fl.zeds[zed_tag] #getting redshift from zed tag

        #loading weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        
        #creating data 
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        MBH += 10. # converting MBH units to M_sol
        
        h = 0.6777 # Hubble parameter

        #converting MBHacc units to M_sol/yr
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr

        x += 10 # converting units of Mstar to M_sol

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
        
        #calculating luminosities with defined function
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
        
        #plotting Lagn vs Mstar at redshift in loop
        plt.scatter(x,lums, alpha = 0.25)#,s=3,c=temps, cmap = cmap,alpha=0.25)
        plt.title('z = ' + str(z))
        plt.ylabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/L_{\odot}yr))$')
        plt.xlabel(r'Stellar Mass ($\rm log_{10}(M*/M_{\odot}))$')
        plt.xlim(8)#,11.5)
        #plt.ylim(-5,15)
        plt.grid()
        fig.savefig(f'figures/LAGN_Mstar/LAGN_Mstar_'+str(int(z))+'.png')
        plt.clf()
        fig.clf()
             
def Lagn_Mstar_grid():
    '''No inputs, creates grid of Lagn vs Mstar at each redshift'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units
 
    flares_dir = '../../../../../../../research/astro/flare' #location of flares
    
    zed_tags = (-1,-2,-3,-4,-5,-6) #tags used to identify redshift in FLARES
    
    #process for loading data and creating arrays is the same as above
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    cmap = mpl.cm.plasma
    
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
        z = fl.zeds[zed_tag]
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        MBH += 10. # convert to M_sol
        h = 0.6777 # Hubble parameter
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr

        x += 10 # units are 1E10 M_sol
        #cutting below mstar = 10**8
        y = MBH
        y_dot = MBHacc
        #cutting below mstar = 10**8
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
        
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
         
        #function is different after this point
        # each redshift plotted into the same figure
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)

        ax.scatter(x, lums, alpha=0.25)
        ax.set_title('z = '+str(z))
        ax.set_ylabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/L_{\odot}))$')
        ax.set_xlabel(r'Stellar Mass ($\rm log_{10}(M*/M_{\odot}))$')
        ax.set_xlim(8,11.5)
        ax.set_ylim(-5,15)
        ax.grid()
        plot_no+=1

    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/LAGN_Mstar/LAGN_Mstar_grid.png')
    plt.clf()
    fig.clf()
