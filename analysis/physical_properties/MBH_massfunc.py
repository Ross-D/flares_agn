def MBH_func():
    '''Creates figure containing six sublots of SMBH massfunction at each redshift'''
    #Importing modules
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import FLARE.plt as fplt
    
    flares_dir = '/research/astro/flare/' # location of flares files
    
    # loading flares
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    
    # setting colour map
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    # creating figure
    fig = plt.figure(figsize=(14,6))
    
    # loading SMBH Mass data
    MBH = fl.load_dataset('BH_Mass', arr_type='Galaxy') # stellar mass of galaxy
    X = MBH

    for i in (-1, -2, -3, -4, -5, -6):
        # geting redshift tag
        tag = fl.tags[i]
        z = fl.zeds[i]
        
        # creating array of data
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        ws, x = np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(X[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(X[halo[ii]][tag]))


        x += 10 # units are 1E10 M_sol
        
        # seting up bins for data
        binw = 0.3
        bins = np.arange(5.35, 9.45, binw)
        b_c = bins[:-1]+binw/2
        
        # binning data (weighted)
        N_weighted, edges = np.histogram(x, bins = bins, weights = ws)

        # calculating mass function
        vol = h = 0.6777
        vol = (4/3)*np.pi*(14/h)**3
        phi = N_weighted/(binw*vol)
        
        bins = np.append(bins, max(bins)+0.5)
        phi = np.append(phi, 1e-100)
        for i in range(len(phi)):
            if phi[i]==0:
                phi[i] = 1e-100
        # plotting at redshift
        plt.plot(bins[:-1]+binw/2, np.log10(phi), c=cmap(norm(z)), label = 'z = '+str(z))

    # Formatting plot
    plt.title('SMBH Mass Function')
    plt.grid()
    plt.ylabel(r'Mass function $\rm log_{10}(\phi/Mpc^{-3}\, dex^{-1})$')
    plt.xlabel(r'SMBH mass $\rm log_{10}(M_{BH}/M_{\odot})$')
    plt.legend()
    plt.ylim(-8, -3)
    plt.xlim(5.5)
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/MBH_massfunc/MBH_single.png')
    plt.clf()
    fig.clf()
    
def Mdot_func():
    '''Creates figure containing six sublots of SMBH accretion rate function at each redshift'''
    # importing modules
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import FLARE.plt as fplt
    import astropy.constants as constants
    import astropy.units as units

    flares_dir = '/research/astro/flare/' # location of flares
    
    # loading in flares
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    
    # setting colour map
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    # creating figure
    fig = plt.figure(figsize=(14,6))
    
    # Loading in data
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black accretion rate
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy') # Black accretion rate

    for i in (-1, -2, -3, -4, -5, -6):
        # redshift tags
        tag = fl.tags[i]
        z = fl.zeds[i]
        
        # reading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        
        # reading in data for redshift
        ws, MBH, MBHdot = np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        
        h = 0.6777 # Hubble parameter
        
        # transforming accretion units
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr
        x = MBHacc
        
        # setting up data bins
        binw = 0.5
        bins = np.arange(-17.75,1.25,binw)
        b_c = bins[:-1]+binw/2
        
        # Binning data
        N_weighted, edges = np.histogram(x, bins = bins, weights = ws)

        #Finding Mass function
        vol = h = 0.6777
        vol = (4/3)*np.pi*(14/h)**3
        phi = N_weighted/(binw*vol)
        
        
        bins = np.append(bins, max(bins)+0.5)
        x = np.array([-18])
        x = np.append(x, bins)
        phi = np.append(phi, 1e-100)
        y = np.zeros(1)
        y = np.append(y, phi)
        phi = y
        bins = x
        for i in range(len(phi)):
            if phi[i]==0:
                phi[i] = 1e-100
        
        # Plotting
        plt.plot(bins[:-1]+binw/2, np.log10(phi), c=cmap(norm(z)), label = 'z = '+str(z))

    # Formatting plot
    plt.title('Accretion Rate "Mass" Function')
    plt.grid()
    plt.ylim(-8, -2)
    plt.xlim(-18, 1)
    plt.ylabel(r'Mass function $\rm log_{10}(\phi/Mpc^{-3}\, dex^{-1})$')
    plt.xlabel(r'Accretion Rate $\rm log_{10}(\dot{M}_{BH}/M_{\odot}yr^{-1})$')
    plt.legend()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'figures/MBH_massfunc/Mdot_single.png')
    plt.clf()
    fig.clf()
