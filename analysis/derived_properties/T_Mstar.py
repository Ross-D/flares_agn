def Temp_Mstar_individual(wanted_zed_tag):
    '''plots T_BB against Mstar, takes single input, must be string:
    'all' gives histograms for all possible redshifts
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
    
    if wanted_zed_tag == 'all':
        zed_tags = (-1,-2,-3,-4,-5,-6)
    else:
        zed_tags = [int(wanted_zed_tag)]
        
    flares_dir = '../../../../../../../research/astro/flare' # location of flares
    
    #importing flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #defining colormap
    cmap = mpl.cm.plasma
    
    #loading in datasets
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
        z = fl.zeds[zed_tag] # convert zed tag into redshift number

        # retrieving weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        #creating arrays for plotting
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        MBH += 10. # converting units of MBH to M_sol
        h = 0.6777 # Hubble parameter
        # converting units of MBHacc to M_sol/yr
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
        
        #calculating temperatures 
        temps = np.array([])
        for t in range(len(x)):
            temp = ((2.24e9)*((10**y_dot[t])**(0.25))*((10**y[t])**(-0.5)))    # Using blue numbers from bluetides 
            temps = np.append(temps, temp)
        
        #plotting individually
        plt.scatter(x, temps/10**6, alpha=0.25)
        plt.xlim(8)
        plt.ylim(0)
        plt.grid()
        plt.title('z = ' + str(z))
        plt.ylabel(r'Temperature (T$_{AGN}$/$10^6$ K)')
        #plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(6,6))
        plt.xlabel(r'Stellar Mass ($\rm log_{10}(M_{\star}/M_{\odot}))$')
        fig.savefig(f'figures/T_Mstar/T_Mstar_'+str(int(z))+'.png')
        plt.clf()
        fig.clf()

def Temp_Mstar_grid():
    '''plots T_BB against Mstar, all plots will be saved as a single figure'''
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
    
    zed_tags = (-1,-2,-3,-4,-5,-6) # tags used to denote redshift in flares
    
    #loading flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #defining colormap
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
        z = fl.zeds[zed_tag] #converting zed tag to redshift number

        #loading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        
        #creating data arrays for plotting
        ws, x, MBH, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        MBH += 10. # converting MBH units to M_sol
        h = 0.6777 # Hubble parameter
        # converting MBHacc units to M_sol/yr
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr

        x += 10 # # converting Mstar units to M_sol
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
        
        #calculating temperatures
        temps = np.array([])
        for t in range(len(x)):
            temp = ((2.24e9)*((10**y_dot[t])**(0.25))*((10**y[t])**(-0.5)))    # Using blue numbers from bluetides 
            temps = np.append(temps, temp)
         
        #plotting into single figure 
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)

        ax.scatter(x, temps/10**6, alpha=0.25)
        ax.set_title('z = '+str(z))
        ax.set_ylabel(r'Temperature (T$_{AGN}$/$10^6$ K)')
        ax.set_xlabel(r'Stellar Mass ($\rm log_{10}(M_{*}/M_{\odot}))$')
        ax.set_ylim(0,1.5)
        ax.set_xlim(8, 11.5)
        ax.grid()
        #ax.ticklabel_format(axis = 'y', style = 'sci', scilimits=(6,6))
        
        temp_max.append(max(temps))
        temp_min.append(min(temps))
        plot_no+=1
    
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f'figures/T_Mstar/T_Mstar_grid_temp.png')
    plt.clf()
    fig.clf()
