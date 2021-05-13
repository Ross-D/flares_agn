def L_AGN(Mdot):
      # This function will calculate the bolometric luminosity of an
    # AGN from an input of accretion rate input of log10(Mdot/(M_sol/yr))
    # the output will be log10(Lagn/(M_sol/yr))
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

def L_AGN_erg(Mdot):
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
    Lagn_CGS = Lagn_SI.to(u.erg/u.s)
    Lagn = np.log10(Lagn_CGS/(u.erg/u.s))
    return Lagn

def Lagn_T_individual(wanted_zed_tag):
    '''Plots T_BB against L_AGN individually, takes single input, must be string:
    'all' gives Lagn vs T for all possible redshifts
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
    
    #sorting which zed tags are needed
    if wanted_zed_tag == 'all':
        zed_tags = (-1,-2,-3,-4,-5,-6)
    else:
        zed_tags = [int(wanted_zed_tag)]
        
    flares_dir = '../../../../../../../research/astro/flare' #location of flares
    
    #importing flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #defining colormap
    cmap = mpl.cm.plasma
    
    #loading required datasets
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
        z = fl.zeds[zed_tag] # converting zedtag to redshift number

        #loading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        # creating data arrays for plotting
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

        x += 10 # changing units of Mstar to M_sol
        
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
        
        #calculating luminosities
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
        #calculating big bump temperatures    
        temps = np.array([])
        for t in range(len(x)):
            temp = ((2.24e9)*((10**y_dot[t])**(0.25))*((10**y[t])**(-0.5))) 
            temps = np.append(temps, temp)
        
        #plotting data individually at each redshift
        plt.scatter(lums, temps/10**6, alpha = 0.25)#,s=3,c=temps, cmap = cmap,alpha=0.25)
        plt.title('z = ' + str(z))
        plt.xlabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/L_{\odot}yr))$')
        plt.ylabel(r'Temperature (T$_{AGN}$/$10^6$ K)')
        plt.ylim(0)
        plt.xlim(0)
        plt.grid()
        #plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(6,6))

        fig.savefig(f'figures/LAGN_T/LAGN_T_'+str(int(z))+'.png')
        plt.clf()
        fig.clf()

def Lagn_T_grid():
    '''Plots T_BB against L_AGN (L_sol units) at all redshifts into a single figure, no inputs required'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units

    flares_dir = '../../../../../../../research/astro/flare' #specify location of flares
    
    zed_tags = (-1,-2,-3,-4,-5,-6) #zed tags used to denote redshift in flares
    
    #importing data from flares
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #specifying colormap
    cmap = mpl.cm.plasma
    
    #loading specific datasets from flares
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
        z = fl.zeds[zed_tag] #converting zed tag to redshift number

        # loading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        
        #creating data arrays to use for plotting        
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
        
        #calculating luminosities
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN(y_dot[l])
            lums = np.append(lums, lum)
            
        #calculating temperatures        
        temps = np.array([])
        for t in range(len(x)):
            temp = ((2.24e9)*((10**y_dot[t])**(0.25))*((10**y[t])**(-0.5)))
            temps = np.append(temps, temp)

        #plotting this redshift
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)

        ax.scatter(lums, temps/(10**6), alpha=0.25)
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/L_{\odot}yr))$')
        ax.set_ylabel(r'Temperature (T$_{AGN}$/$10^6$ K)')
        ax.set_ylim(0, 1.5)
        ax.set_xlim(0, 15)
        ax.grid()
        #ax.ticklabel_format(axis = 'y', style = 'sci', scilimits=(6,6))
        plot_no+=1
  
    plt.tight_layout(pad=0.5, w_pad=1.25, h_pad=0.5)
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/LAGN_T/LAGN_T_grid.png')
    plt.clf()
    fig.clf()
    
def Lagn_T_grid_erg():
    '''Plots T_BB against L_AGN (erg/s units) at all redshifts into a single figure, no inputs required, prcoess is the same as in above function'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units

    flares_dir = '../../../../../../../research/astro/flare'
    
    zed_tags = (-1,-2,-3,-4,-5,-6)
    
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
        
        #only diffference from above function, luminosities calculated in units of erg/s
        lums = np.array([])
        for l in range(len(x)):
            lum = L_AGN_erg(y_dot[l])
            lums = np.append(lums, lum)
            
        temps = np.array([])
        for t in range(len(x)):
            temp = ((2.24e9)*((10**y_dot[t])**(0.25))*((10**y[t])**(-0.5)))    # Using blue numbers from bluetides 
            temps = np.append(temps, temp)

        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)

        ax.scatter(lums, temps/(10**6), alpha=0.25)
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r'AGN Luminosity ($\rm log_{10}(L_{AGN}/ergs^{-1}))$')
        ax.set_ylabel(r'Temperature (T$_{AGN}$/$10^6$ K)')
        ax.set_ylim(0, 1.5)
        ax.set_xlim(32.5, 48)
        ax.grid()
        #ax.ticklabel_format(axis = 'y', style = 'sci', scilimits=(6,6))
        plot_no+=1

    plt.tight_layout(pad=0.5, w_pad=1.25, h_pad=0.5)
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/LAGN_T/LAGN_T_grid_erg.png')
    plt.clf()
    fig.clf()
