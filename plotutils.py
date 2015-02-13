import re, os, sys

def scatter1(filekeyword, xvalues, yvalues, xmin=None, ymin=None, xmax=None, ymax=None, logy=False, logx=False, xlab="", ylab="", format="pdf", colors=["blue", "red", "purple"]):
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from matplotlib.gridspec import GridSpec
    import pylab as P
    


    """If the user didn't specify xmin and ymin, then let's
        auto-determine these values from the data."""
    if xmin == None:
        xmin = min(xvalues[0])
    if ymin == None:
        ymin = min(yvalues[0])
    if xmax == None:
        xmax = max(xvalues[0])
    if ymax == None:
        ymax = 1.05 * max(yvalues[0])
    
    fig1 = plt.figure(figsize=[8,6])
    gs = GridSpec(100, 100, bottom=0.10,left=0.15,right=0.85)
    ax1 = fig1.add_subplot(gs[15:100,0:85])
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    if logx:
        ax1.set_xscale("log")
    if logy:
        ax1.set_yscale("log")
    #ax1.set_yscale("log")
    
    ax1.scatter(xvalues[0], yvalues[0], color='blue', alpha=0.3, edgecolors='none', s=5)
    ax1.set_xlabel(xlab, fontsize=14)
    ax1.set_ylabel(ylab, fontsize=14)
    
    """Also plot the other data series"""
    if xvalues.__len__() > 1:
        for ii in range(1, xvalues.__len__()):
            ax1.scatter(xvalues[ii], yvalues[ii], color=colors[ii], edgecolors='none', s=5)
    
    """Set ticks and axes"""
    ax1.set_xticks([xmin,xmax])
    ax1.set_yticks([ymin,ymax])
    ax1.set_xticklabels([xmin.__str__(), xmax.__str__()])
    ax1.set_yticklabels([ymin.__str__(), ymax.__str__()] )
    ax1.tick_params(axis='x',length=3,width=1.5)
    
    """Top histogram"""
    ax2 = fig1.add_subplot(gs[0:14,0:85])
    ax2.set_xlim(xmin,xmax)
    if logx:
        n, bins, patches = ax2.hist(xvalues[0], bins = 10 ** np.linspace(np.log10(xmin), np.log10(xmax), 50))
        ax2.set_xscale("log")
    else:
        n, bins, patches = ax2.hist(xvalues[0], bins = np.linspace(xmin, xmax, 50) )
    
    """Right-side histogram"""    
    ax3 = fig1.add_subplot(gs[15:100,86:100])
    ax3.set_ylim(ymin,ymax)
    if logy:
        n, bins, patches = ax3.hist(yvalues[0], orientation="horizontal", bins = 10 ** np.linspace(np.log10(ymin), np.log10(ymax), 50))
        ax3.set_yscale("log")
    else:
        n, bins, patches = ax3.hist(yvalues[0], orientation="horizontal", bins = np.linspace(ymin, ymax, 50) )

    """Labels for the histograms"""
    ax2.set_frame_on(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)

    ax3.set_frame_on(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)


    """Finally, save the image."""
    fig1.savefig(filekeyword + "." + format, format=format, dpi=600)


