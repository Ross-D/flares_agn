def MBH_hist_single_line():
    '''Plots histogram of SMBH masses at all redshifts onto a single plot'''    
    #importing modules
    import numpy as np
    import pandas as pd
    import numpy as np
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    
    flares_dir = '../../../../../../../research/astro/flare' #location of flares directory
    zed_tags = (-1,-2,-3,-4,-5,-6) #tags used to specify redshifts in flares
    
    #importing flares data set
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    #creating and normalising colour map
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    #extracting stellar mass and black hole mass for the required sim and redshift
    Mstar = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBH = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy
    
    X = Mstar
    Y = MBH

    fig = plt.figure(figsize=(14,6))
    #plot_no = 1
    
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag]
        #creating data arrays
        x, y = np.array([]), np.array([])
        for ii in range(len(halo)):
            y = np.append(y, np.log10(Y[halo[ii]][tag]))
            x = np.append(x, np.log10(X[halo[ii]][tag]))
        #converting MBH and Mstar to units of M_sol
        x += 10 
        y += 10 
        #cuttting out unseeded galaxies
        y2 = []
        for i in range(len(y)):
            if y[i] == -np.inf:
                continue
            else:
                y2.append(y[i])
        
        #binning SMBHs into mass bins 
        bins = np.arange(5, 10.15, 0.3)
        y_binned = np.histogram(y2, bins = bins)
        bincen = []
        for i in range(len(y_binned[1])-1):
            cen = y_binned[1][i] + (y_binned[1][i+1] - y_binned[1][i])/2
            bincen.append(cen)
        #plotting histogram
        plt.plot(bincen, y_binned[0], c = cmap(norm(z)), label = 'z = '+str(int(z)))#[0], y_binned[1])

    plt.ylabel('Frequency')
    plt.yscale('log')
    plt.xlabel(r'$\rm log_{10}(M_{BH}/M_{\odot})$')
    plt.legend()
    plt.xlim(5.5, 9.5)
    plt.grid()
    plt.title('Total Simulated SMBH by Mass')
    fig.savefig(f'figures/MBH_hist/hist_single_line.png')
    plt.clf()
    fig.clf()
