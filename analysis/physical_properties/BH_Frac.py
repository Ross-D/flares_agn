def plotBH_frac_single():
    '''Plots the fraction of galaxies that contain a SMBH binned into masses at each redshift onto a single plot'''
    #importing modules
    import numpy as np
    import pandas as pd
    import numpy as np
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    
    flares_dir = '../../../../../../../research/astro/flare' # location of FLARES directory
    zed_tags = (-1,-2,-3,-4,-5,-6) # need all redshifts
    
    # loading in flares
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    
    # setting up colourmap
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    # loading data from flares
    Mstar = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBH = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy
    
    X = Mstar
    Y = MBH
    
    # creating figure
    fig = plt.figure(figsize=(14,6))
    
    # plotting for each redshift
    for zed_tag in zed_tags:
        
        # getting redshift tag
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag]
        
        # importing data
        x, y = np.array([]), np.array([])
        for ii in range(len(halo)):
            y = np.append(y, np.log10(Y[halo[ii]][tag]))
            x = np.append(x, np.log10(X[halo[ii]][tag]))
            
        x += 10 # units are 1E10 M_sol
        y += 10 # units are 1E10 M_sol
        
        # creating set of just galaxies with AGN
        x_with = []
        for i in range(len(y)):
            if y[i] == -np.inf:
                continue
            else:
                x_with.append(x[i])
        
        # binning data
        bins = 18 # int((max(x)-min(x))*4)+1
        gals = np.histogram(x, bins = bins)
        gals_with = np.histogram(x_with, bins = bins)
        
        # fraction of galaxies with SMBH in each bin
        fracs = gals_with[0]/gals[0]
        
        # plotting data for this redshift
        linewidth = 0.8*(gals[1][1]-gals[1][0])
        bincen = []
        for i in range(len(fracs)):
            cen = (gals[1][i]+gals[1][i+1])/2
            bincen.append(cen)
        plt.plot(bincen, fracs, c=cmap(norm(z)), label = 'z = '+str(int(z)))
    
    # formatting plots
    plt.legend()
    plt.ylabel(r'Fraction with SMBH')
    plt.xlabel(r'Stellar Mass ($\rm log_{10}(M_{\star}/M_{\odot})$)')
    plt.xlim(7.2, 11.5)
    #plt.ylim(0.8, 1.01)
    plt.grid()
    plt.savefig(f'figures/BH_Frac/single.png')
    plt.clf()
