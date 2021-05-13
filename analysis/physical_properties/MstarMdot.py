def Mstar_Mdot_grid_MBH():
    '''plots Mstar against Mdot at all redshifts in a single figure'''
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.units as units
    import astropy.constants as constants

    flares_dir = '../../../../../../../research/astro/flare' #location of flares
    
    zed_tags = (-1,-2,-3,-4,-5,-6) # tags used to specify redshift in flares
    
    # importing flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #specifying and normalising colormap
    cmap = mpl.cm.cividis_r
    norm = mpl.colors.Normalize(vmin=5, vmax=9.5)
    
    #importing datasets from flares
    Mstar_ = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy')

    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag] # converting zed tag to redshift number
        #loading in weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        #creating data arrays        
        ws, x, y, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(Mstar_[halo[ii]][tag]))
            y = np.append(y, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))

        x += 10 # converting units of Mstar to M_sol
        y += 10 # converting units of MBH to M_sol
        
        h = 0.6777 # Hubble parameter
        
        #converting units of MBHacc to M_sol/yr
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr
        y_dot = MBHacc
        

        # cutting out galaxies below 10^8 Msol
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
         # removing unseeded galaxies
        xcut, ycut, y_dotcut = np.array([]), np.array([]), np.array([])
        for ii in range(len(y)):
            if y[ii] == -np.inf:
                continue
            else:
                xcut = np.append(xcut, x[ii])
                ycut = np.append(ycut, y[ii])
                y_dotcut = np.append(y_dotcut, y_dot[ii])
        x = xcut
        y = ycut
        y_dot = y_dotcut
        
        #plotting MBHacc against Mstar, colouring by MBH
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(x, y_dot, s=3,c=y, cmap = cmap, norm = norm, alpha=0.5)
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r'Stellar Mass ($\rm log_{10}(M_{\star}/M_{\odot}))$')
        ax.set_ylabel(r'Accretion Rate, $log_{10}(\dot{M}_{BH}/M_\odot yr^{-1})$')
        ax.set_xlim(8, 11.5)
        ax.set_ylim(-18, 1)
        ax.grid()
        plot_no+=1

    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.subplots_adjust(right=0.85)
    # adding colourbar
    cbar_ax = fig.add_axes([0.875, 0.1, 0.01, 0.8])
    cbar = fig.colorbar(cm.ScalarMappable(cmap = cmap, norm = norm), cax = cbar_ax)
    cbar.set_label(r'SMBH Mass ($\rm log_{10}(M_{BH}/M_{\odot}))$', rotation=270, fontsize = 7, labelpad = 8)
    cbar.ax.tick_params(labelsize=7)

    fig.savefig(f'figures/MstarMdot/MstarMdot_grid_MBH.png')
    plt.clf()
    fig.clf()
    
def Mstar_Mdot_individual_MBH(wanted_zed_tag):
    '''plots Mstar against Mdot at all redshifts individually,
     takes single input, must be string:
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
    import astropy.units as units
    import astropy.constants as constants

    flares_dir = '../../../../../../../research/astro/flare'#location of flares
    
    #sorting needed zed tags
    if wanted_zed_tag == 'all':
        zed_tags = (-1,-2,-3,-4,-5,-6)
    else:
        zed_tags = [int(wanted_zed_tag)]
    
    #importing in flares data
    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos

    #defining and normalising colourmap
    cmap = mpl.cm.cividis_r
    norm = mpl.colors.Normalize(vmin = 5, vmax = 9.5)
    
    #loading in datasets
    Mstar_ = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy')

    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag] # converting zed tag to redshift number

        #loading weights
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        
        #creating data arrays
        ws, x, y, MBHdot = np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(Mstar_[halo[ii]][tag]))
            y = np.append(y, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
        
        #converting units of MBH and Mstar to M_sol
        x += 10 
        y += 10 
        
        h = 0.6777 # Hubble parameter
        
        # converting units of MBHacc to M_sol/yr
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr
        y_dot = MBHacc
        
        # cutting out galaxies below Mstar = 10^8
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
        
        #removing unseeded galaxies
        xcut, ycut, y_dotcut = np.array([]), np.array([]), np.array([])
        for ii in range(len(y)):
            if y[ii] == -np.inf:
                continue
            else:
                xcut = np.append(xcut, x[ii])
                ycut = np.append(ycut, y[ii])
                y_dotcut = np.append(y_dotcut, y_dot[ii])
        x = xcut
        y = ycut
        y_dot = y_dotcut

        #plotting mbhacc against mstar, colouring by MBH
        plt.scatter(x, y_dot,s=3, c=y, cmap = cmap, norm = norm, alpha=0.9)
        plt.title('z = ' + str(z))
        plt.ylabel(r'Accretion Rate ($\rm log_{10}(\dot{M}/M_{\odot}yr))$')
        plt.xlabel(r'Stellar Mass ($log_{10}(M_\star/M_\odot$))')
        plt.xlim(8, 11.5)
        plt.ylim(-18, 1)
        plt.grid()
        plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
        fig.subplots_adjust(right=0.85)
        cbar = fig.colorbar(cm.ScalarMappable(cmap = cmap, norm = norm))
        #cbar.set_ticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.3])
        cbar.set_label(r'SMBH Mass ($\rm log_{10}(M_{BH}/M_{\odot}))$', rotation=270, fontsize = 7, labelpad = 8.5)
        cbar.ax.tick_params(labelsize=7)
        fig.savefig(f'figures/MstarMdot/MstarMdot_'+str(int(z))+'_MBH.png')
        plt.clf()
        fig.clf()
