def L_AGN(Mdot):
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

def Temp(m, mdot):
    return ((2.24e9)*((10**mdot)**(0.25))*((10**m)**(-0.5)))

def L_FUV_fit(l_bol, temp):
    import numpy as np
    x = temp/1e5
    a = -1.768e-3
    b = 4.901e-1
    c = -3.957e-3
    d = 6.158e-2
    e = 2.123
    f = 5.954e-2
    g = 2.547
    ratio = (g*((a*(x**b))+(c/(d+(x**e)))+f))
    if ratio > 0:
        Ratio = ratio
    else:
        Ratio = 1e-8
    return np.log10(Ratio*(10**l_bol))

def retrieve_data(): 
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import astropy.constants as constants
    import astropy.units as units


    zed_tags = (-1,-2,-3,-4,-5,-6)
    #zed_tags = [-1]        

    flares_dir = '../../../../../../../research/astro/flare'

    fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
    halo = fl.halos
    cmap = mpl.cm.plasma

    Mstar_ = fl.load_dataset('Mstar_30', arr_type='Galaxy') # stellar mass of galaxy
    MBHdot_ = fl.load_dataset('BH_Mdot', arr_type='Galaxy') # Black accretion rate
    MBH_ = fl.load_dataset('BH_Mass', arr_type='Galaxy') # Black hole mass of galaxy
    Lum_Gal = fl.load_dataset("BPASS_2.2.1/Chabrier300/Luminosity/Intrinsic/FUV", arr_type='Galaxy') # Black hole mass of galaxy
    Data = {}
    for zed_tag in zed_tags:
        tag = fl.tags[zed_tag]
        z = fl.zeds[zed_tag]
        df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        ws, Mstar, MBH, MBHdot, Gal_FUV = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(MBH_[halo[ii]][tag]))*weights[ii])
            Mstar = np.append(Mstar, np.log10(Mstar_[halo[ii]][tag]))
            MBH = np.append(MBH, np.log10(MBH_[halo[ii]][tag]))
            MBHdot = np.append(MBHdot, np.log10(MBHdot_[halo[ii]][tag]))
            Gal_FUV = np.append(Gal_FUV, np.log10(Lum_Gal[halo[ii]][tag]))
        Mstar += 10.  # convert to M_sol
        MBH += 10. # convert to M_sol
        h = 0.6777 # Hubble parameter
        
        MBHacc = MBHdot + np.log10(h*6.445909132449984E23) # g/s
        MBHacc -= np.log10(constants.M_sun.to('g').value) # convert to M_sol/s
        MBHacc += np.log10(units.yr.to('s')) # convert to M_sol/yr
        MBH_cut, MBHacc_cut, Gal_FUV_cut, Mstar_cut, ws_cut = [], [], [], [], []
        for i in range(len(MBH)):
            if MBH[i] == -np.inf:
                continue
            elif MBHacc[i] == -np.inf:
                continue
            else:
                MBH_cut.append(MBH[i])
                MBHacc_cut.append(MBHacc[i])
                Gal_FUV_cut.append(Gal_FUV[i])
                Mstar_cut.append(Mstar[i])
                ws_cut.append(ws[i])
        
        MBH = np.asarray(MBH_cut)
        MBHacc = np.asarray(MBHacc_cut)
        Gal_FUV = np.asarray(Gal_FUV_cut)
        Mstar = np.asarray(Mstar_cut)
        ws = np.asarray(ws_cut)
        
        AGN_Lbol, AGN_FUV, AGN_T = np.zeros(len(MBHacc)), np.zeros(len(MBHacc)), np.zeros(len(MBHacc))
        for iii in range(len(MBHacc)):
            L_agn = L_AGN(MBHacc[iii])
            T_agn = Temp(MBH[iii], MBHacc[iii])
            L_agn_FUV = L_FUV_fit(L_agn, T_agn)
            AGN_Lbol[iii] = L_agn + 18.27875360095283
            AGN_FUV[iii] = L_agn_FUV + 18.27875360095283
            AGN_T[iii] = T_agn
            
        data_at_z = {}
        data_at_z.update({"AGN_Lbol":AGN_Lbol})
        data_at_z.update({"AGN_FUV":AGN_FUV})
        data_at_z.update({"AGN_T":AGN_T})
        data_at_z.update({"Gal_FUV":Gal_FUV})
        data_at_z.update({"MBH":MBH})
        data_at_z.update({"MBHacc":MBHacc})
        data_at_z.update({"Mstar":Mstar})
        data_at_z.update({"Weights":ws})
        Data.update({str(z):data_at_z})
    return Data
