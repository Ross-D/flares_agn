def Lagn_FUV_Lgal_FUV():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        linecolor = "#00008b"
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, AGN_FUV,ws,bins,quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
        
        x=np.linspace(0,35)
        plt.figure(figsize=(14, 8))
        plt.scatter(Gal_FUV, AGN_FUV, alpha = 0.25)
        plt.plot(x, x, "k--", label = r"L$_{AGN,FUV}$ = L$_{Gal,FUV}$")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.legend(loc = "lower right")
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        plt.ylim(min(AGN_FUV)-0.5, max(AGN_FUV)+0.5)
        plt.xlim(28, max(Gal_FUV)+0.5)  
        plt.grid()
        plt.savefig(f"Plots/Lfuv_AGN_GAL/Lfuv_AGN_GAL_"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        linecolor = "#00008b"
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, AGN_FUV,ws,bins,quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
        
        x=np.linspace(0,35)
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Gal_FUV, AGN_FUV, alpha = 0.25)
        ax.plot(x, x, "k--", label = r"L$_{AGN,FUV}$ = L$_{Gal,FUV}$")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.legend(loc = "lower right")
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.grid()
        ax.set_ylim(8, 31.6)
        ax.set_xlim(28, 31.6)  
        plot_no+=1
    x=np.linspace(0,35)  
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lfuv_AGN_GAL/Lfuv_AGN_GAL_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        x=np.linspace(0,35)
        
        ws = Data[z]["Weights"]
        linecolor = "#00008b"
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, AGN_FUV,ws,bins,quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
             
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Gal_FUV, AGN_FUV, alpha = 0.25)
        ax.plot(x, x, "k--", label = r"L$_{AGN,FUV}$ = L$_{Gal,FUV}$")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.legend(loc = "lower right")
        ax.grid()
        ax.set_ylim(28, 31.6)
        ax.set_xlim(28, 31.2)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lfuv_AGN_GAL/Lfuv_AGN_GAL_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def FUV_ratio_Lgal():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    
    linecolor = "#00008b"
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        x = np.linspace(min(Gal_FUV)-0.5, max(Gal_FUV)+0.5)
        y = np.zeros(len(x))
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
        
        plt.figure(figsize=(14, 8))
        plt.scatter(Gal_FUV, Ratio, alpha = 0.25)
        plt.plot(x, y, "k--", label = r"L$_{AGN,FUV}$ = L$_{Gal,FUV}$")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        #plt.ylim(-14.2, 2)
        plt.xlim(28, max(x))  
        plt.grid()
        
        plt.legend(loc = "lower right")
        
        plt.savefig(f"Plots/FUV_ratio_Lgal/FUV_ratio_Lgal_"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(25.1, 31.2)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Gal_FUV, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "lower right")
        
        ax.set_ylim(-20, 2)
        ax.set_xlim(28, 31.2)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_Lgal/FUV_ratio_Lgal_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(25.1, 31.2)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(27.9, 31.3, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Gal_FUV, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Gal_FUV, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Gal_FUV, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(L$_{Gal,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "upper left")
        
        ax.set_ylim(-1, 2)
        ax.set_xlim(28, 31.2)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_Lgal/FUV_ratio_Lgal_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def AGN_Lgal_hist():
    import numpy as np
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    fig = plt.figure(figsize=(14,6))
    for z in zeds:
        BH_Mass = Data[z]["MBH"]
        Gal_FUV = Data[z]["Gal_FUV"]
        Gal_FUV_cut = np.array([])

        for i in range(len(BH_Mass)):
            if BH_Mass[i] != np.inf:
                Gal_FUV_cut = np.append(Gal_FUV_cut, Gal_FUV[i])
            else:
                continue
        
        bins = np.arange(25, 31.5, 0.2)
        Gal_FUV_binned = np.histogram(Gal_FUV_cut, bins = bins)
        bincen = []
        for i in range(len(Gal_FUV_binned[1])-1):
            cen = Gal_FUV_binned[1][i] + (Gal_FUV_binned[1][i+1] - Gal_FUV_binned[1][i])/2
            bincen.append(cen)

        bincen = np.asarray(bincen)
        plt.plot(bincen, Gal_FUV_binned[0], c = cmap(norm(float(z))), label = 'z = '+str(int(float(z))))#[0], y_binned[1])
    plt.ylabel('Number of AGN')
    plt.yscale('log')
    plt.xlabel(r'log$_{10}$(L$_{GAL,FUV}$/erg s$^{-1}$ Hz^${-1}$)')
    plt.xlim(25, 31.2)
    plt.legend()
    plt.grid()
    plt.title(r'Number of AGN as a function of L$_{GAL,FUV}$')
    plt.savefig(f'Plots/Lgal_AGN_Hist/Lgal_AGN_Hist.png')
    plt.clf()
    fig.clf() 
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def Lagn_FUV_Mstar():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    linecolor = "#ff6700"
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        Mstar = Data[z]["Mstar"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
                
        plt.figure(figsize=(14, 8))
        plt.scatter(Mstar, AGN_FUV, alpha = 0.25)
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        #plt.ylim(14.2, 31.6)
        plt.xlim(8)  
        plt.grid()
        
        plt.legend(loc = "lower right")
        
        plt.savefig(f"Plots/Lagn_FUV_Mstar/Lagn_FUV_Mstar"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Mstar = Data[z]["Mstar"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Mstar, AGN_FUV, alpha = 0.25)
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.grid()
        
        ax.legend(loc = "lower right")
        
        ax.set_ylim(10, 32)
        ax.set_xlim(8, 12)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lagn_FUV_Mstar/Lagn_FUV_Mstar_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Mstar = Data[z]["Mstar"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Mstar, AGN_FUV, alpha = 0.25)
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.grid()
        
        plt.legend(loc = "lower right")
        
        ax.set_ylim(28, 32)
        ax.set_xlim(8, 12)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lagn_FUV_Mstar/Lagn_FUV_Mstar_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def FUV_ratio_Mstar():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    
    linecolor = "#39FF14"
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        Mstar = Data[z]["Mstar"]
        x = np.linspace(min(Mstar)-0.5, max(Mstar)+0.5)
        y = np.zeros(len(x))
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
            
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
        
        plt.figure(figsize=(14, 8))
        plt.scatter(Mstar, Ratio, alpha = 0.25)
        plt.plot(x, y, "k--")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        #plt.ylim(14.2, 31.6)
        plt.xlim(8, max(x))
        
        plt.legend(loc = "lower right")
        
        plt.grid()
        plt.savefig(f"Plots/FUV_ratio_Mstar/FUV_ratio_Mstar_"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        Mstar = Data[z]["Mstar"]
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(7, 12)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
                
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Mstar, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "lower right")
        
        ax.set_ylim(-17, 2)
        ax.set_xlim(8, 11.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_Mstar/FUV_ratio_Mstar_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        Mstar = Data[z]["Mstar"]
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(7, 12)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(7.9, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(Mstar, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(Mstar, bins=bins)
        Ns = N>10
                
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(Mstar, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{\star}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "upper right")
        
        ax.set_ylim(-1, 2)
        ax.set_xlim(8, 11.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_Mstar/FUV_ratio_Mstar_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def Lagn_FUV_MBH():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    
    linecolor = "#FF0000"
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        MBH = Data[z]["MBH"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        plt.figure(figsize=(14, 8))
        plt.scatter(MBH, AGN_FUV, alpha = 0.25)
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        #plt.ylim(14.2, 31.6)
        plt.xlim(5.16)  
        plt.grid()
        
        plt.legend(loc = "lower right")
        
        plt.savefig(f"Plots/Lagn_FUV_MBH/Lagn_FUV_MBH_"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        MBH = Data[z]["MBH"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(MBH, AGN_FUV, alpha = 0.25)
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.grid()
        
        ax.legend(loc = "lower right")
        
        ax.set_ylim(5, 32)
        ax.set_xlim(5.16, 9.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lagn_FUV_MBH/Lagn_FUV_MBH_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        MBH = Data[z]["MBH"]
        AGN_FUV = Data[z]["AGN_FUV"]
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, AGN_FUV, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(MBH, AGN_FUV, alpha = 0.25)
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/erg s$^{-1}$ Hz$^{-1}$)")
        ax.grid()
        
        ax.legend(loc = "upper left")
        
        ax.set_ylim(28, 32)
        ax.set_xlim(5.16, 9.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/Lagn_FUV_MBH/Lagn_FUV_MBH_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
def FUV_ratio_MBH():
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import flares
    
    linecolor = "#800080"
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        MBH = Data[z]["MBH"]
        x = np.linspace(5, max(MBH)+0.5)
        y = np.zeros(len(x))
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
            
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        plt.figure(figsize=(14, 8))
        plt.scatter(MBH, Ratio, alpha = 0.25)
        plt.plot(x, y, "k--")
        
        plt.plot(bincen, out[:,1], c=linecolor, ls= ':')
        plt.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        plt.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        plt.title("z="+str(z))
        plt.xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        plt.ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        #plt.ylim(14.2, 31.6)
        plt.xlim(5.16, max(x))  
        plt.grid()
        
        plt.legend(loc = "lower right")
        
        plt.savefig(f"Plots/FUV_ratio_MBH/FUV_ratio_MBH_"+str(int(float(z)))+".png")
        plt.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        MBH = Data[z]["MBH"]
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(5, 12)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(MBH, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "lower right")
        
        ax.set_ylim(-17, 2)
        ax.set_xlim(5.16, 9.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_MBH/FUV_ratio_MBH_grid.png")
    plt.clf()
    fig.clf()
    
    fig = plt.figure(figsize=(14,6))
    plot_no = 1
    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        MBH = Data[z]["MBH"]
        Ratio = np.zeros(len(Gal_FUV))
        for i in range(len(Ratio)):
            Ratio[i] = AGN_FUV[i]-Gal_FUV[i]
        
        x = np.linspace(5, 12)
        y = np.zeros(len(x))
        
        ws = Data[z]["Weights"]
        quantiles = [0.84,0.50,0.16] # quantiles for range
        bins = np.arange(4.8, 12, 0.2) # x-coordinate bins, in this case stellar mass
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(MBH, Ratio, ws, bins, quantiles)
        N, bin_edges = np.histogram(MBH, bins=bins)
        Ns = N>10
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.scatter(MBH, Ratio, alpha = 0.25)
        ax.plot(x, y, "k--")
        
        ax.plot(bincen, out[:,1], c=linecolor, ls= ':')
        ax.plot(bincen[Ns], out[:,1][Ns], c=linecolor, label = r'Weighted median')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=linecolor, alpha = 0.2)
        
        ax.set_title('z = '+str(z))
        ax.set_xlabel(r"log$_{10}$(M$_{BH}$/M$_\odot$)")
        ax.set_ylabel(r"log$_{10}$(L$_{AGN,FUV}$/L$_{Gal,FUV}$)")
        ax.grid()
        
        ax.legend(loc = "upper left")
        
        ax.set_ylim(-1, 2)
        ax.set_xlim(5.16, 9.5)  
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_ratio_MBH/FUV_ratio_MBH_grid_zoom.png")
    plt.clf()
    fig.clf()
#-------------------------------------------------------------------------------------------------------------
def Lfunc_single():
    #Importing modules
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import FLARE.plt as fplt
    
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    
    # setting colour map
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    # AGN L_FUNC
    # creating figure
    fig = plt.figure(figsize=(14,6))

    for z in zeds:
        AGN_FUV = Data[z]["AGN_FUV"]
        ws = Data[z]["Weights"]
        # seting up bins for data
        binw = 0.4
        bins = np.arange(8, 32, binw)
        b_c = bins[:-1]+binw/2
        
        # binning data (weighted)
        N_weighted, edges = np.histogram(AGN_FUV, bins = bins, weights = ws)

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
        plt.plot(bins[:-1]+binw/2, np.log10(phi), c=cmap(norm(float(z))), label = 'z = '+str(float(z)))

    # Formatting plot
    plt.title(r'AGN FUV Luminosity Function')
    plt.grid()
    plt.ylabel(r'Luminosity function $\rm log_{10}(\phi_{AGN}/Mpc^{-3}\, dex^{-1})$')
    plt.xlabel(r' $\rm log_{10}(L_{AGN,FUV}/erg s^{-1} Hz^{-1})$')
    plt.ylim(-8, -1)
    plt.xlim(28, 31.5)
    plt.legend()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'Plots/FUV_Lfunc/AGN_FUV_Lfunc.png')
    plt.clf()
    fig.clf()
    
    
    # GALAXY L_FUNC
    # creating figure
    fig = plt.figure(figsize=(14,6))

    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        ws = Data[z]["Weights"]
        # seting up bins for data
        binw = 0.4
        bins = np.arange(8, 32, binw)
        b_c = bins[:-1]+binw/2
        
        # binning data (weighted)
        N_weighted, edges = np.histogram(Gal_FUV, bins = bins, weights = ws)

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
        plt.plot(bins[:-1]+binw/2, np.log10(phi), c=cmap(norm(float(z))), label = 'z = '+str(float(z)))

    # Formatting plot
    plt.title(r'Galaxy FUV Luminosity Function')
    plt.grid()
    plt.ylabel(r'Luminosity function $\rm log_{10}(\phi_{Gal}/Mpc^{-3}\, dex^{-1})$')
    plt.xlabel(r' $\rm log_{10}(L_{Gal,FUV}/erg s^{-1} Hz^{-1})$')
    plt.ylim(-8, -1)
    plt.xlim(28, 31.5)
    plt.legend()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'Plots/FUV_Lfunc/Gal_FUV_Lfunc.png')
    plt.clf()
    fig.clf()
    
    # Total L_FUNC
    # creating figure
    fig = plt.figure(figsize=(14,6))

    for z in zeds:
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        Total_FUV = np.zeros(len(Gal_FUV))
        for i in range(len(Gal_FUV)):
            Total_FUV[i] = np.log10((10**Gal_FUV[i])+(10**AGN_FUV[i]))
        ws = Data[z]["Weights"]
        # seting up bins for data
        binw = 0.4
        bins = np.arange(8, 32, binw)
        b_c = bins[:-1]+binw/2
        
        # binning data (weighted)
        N_weighted, edges = np.histogram(Total_FUV, bins = bins, weights = ws)

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
        plt.plot(bins[:-1]+binw/2, np.log10(phi), c=cmap(norm(float(z))), label = 'z = '+str(float(z)))

    # Formatting plot
    plt.title(r'Total FUV Luminosity Function')
    plt.grid()
    plt.ylabel(r'Luminosity function $\rm log_{10}(\phi_{Total}/Mpc^{-3}\, dex^{-1})$')
    plt.xlabel(r' $\rm log_{10}(L_{Gal+AGN,FUV}/erg s^{-1} Hz^{-1})$')
    plt.ylim(-8, -1)
    plt.xlim(28, 31.5)  
    plt.legend()
    fig.subplots_adjust(right=0.85)
    fig.savefig(f'Plots/FUV_Lfunc/Total_FUV_Lfunc.png')
    plt.clf()
    fig.clf()
    
def Lfunc_comparison():
    #Importing modules
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.cm as cm
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import flares
    import FLARE.plt as fplt
    zeds = ("5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
    
    # setting colour map
    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=5., vmax=10.)
    
    # AGN L_FUNC
    # creating figure
    fig = plt.figure(figsize=(14,6))
    plot_no = 1

    for z in zeds:
        
        
        Gal_FUV = Data[z]["Gal_FUV"]
        AGN_FUV = Data[z]["AGN_FUV"]
        Total_FUV = np.zeros(len(Gal_FUV))
        for i in range(len(Gal_FUV)):
            Total_FUV[i] = np.log10((10**Gal_FUV[i])+(10**AGN_FUV[i]))
        ws = Data[z]["Weights"]
        # seting up bins for data
        binw = 0.4
        bins = np.arange(8, 32, binw)
        b_c = bins[:-1]+binw/2
        
        # binning data (weighted)
        AGN_weighted, edges = np.histogram(AGN_FUV, bins = bins, weights = ws)
        Gal_weighted, edges = np.histogram(Gal_FUV, bins = bins, weights = ws)
        Total_weighted, edges = np.histogram(Total_FUV, bins = bins, weights = ws)
        # calculating mass function
        vol = h = 0.6777
        vol = (4/3)*np.pi*(14/h)**3
        phi_AGN = AGN_weighted/(binw*vol)
        phi_Gal = Gal_weighted/(binw*vol)
        phi_Total = Total_weighted/(binw*vol)
        bins = np.append(bins, max(bins)+0.5)
        
        phi_AGN = np.append(phi_AGN, 1e-100)
        phi_Gal = np.append(phi_Gal, 1e-100)
        phi_Total = np.append(phi_Total, 1e-100)
        
        for i in range(len(phi_AGN)):
            if phi_AGN[i]==0:
                phi_AGN[i] = 1e-100
            if phi_Total[i]==0:
                phi_Total[i] = 1e-100
            if phi_Gal[i]==0:
                phi_Gal[i] = 1e-100
            else:
                continue
        
        figloc = '23'+str(plot_no)
        ax = fig.add_subplot(figloc)
        ax.plot(bins[:-1]+binw/2, np.log10(phi_Gal), label = "Galaxy")
        ax.plot(bins[:-1]+binw/2, np.log10(phi_AGN), label = "AGN")
        ax.plot(bins[:-1]+binw/2, np.log10(phi_Total), label = "Total")
        ax.set_title('z = '+str(z))
        ax.set_ylabel(r"Luminosity function $\rm log_{10}(\phi/Mpc^{-3}\, dex^{-1})$")
        ax.set_xlabel(r' $\rm log_{10}(L_{FUV}/erg s^{-1} Hz^{-1})$')
        ax.grid()
        ax.legend()
        ax.set_ylim(-8, -1)
        ax.set_xlim(28, 31.5)   
        plot_no+=1
    plt.tight_layout(pad=0.5, w_pad=1, h_pad=0.5)
    fig.savefig(f"Plots/FUV_Lfunc/Compare_FUV_Lfunc_grid.png")
    plt.clf()
    fig.clf()
