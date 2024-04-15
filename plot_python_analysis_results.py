import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
from matplotlib.colors import *
from CLASSES import *
import math


'''
##########################################
##                                                                                         ##
##                                      Plot Ticks                                     ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to plot_python_analysis_determine_ticks     ##########

def plot_python_analysis_determine_ticks(min, max, number_of_ticks=5, tickformat="%g", scale="linear"):

	range						= max - min
	log_range				= 0

	# handle log scale
	if "log" in scale:
		min_oom			= math.floor(math.log10(min))
		max_oom			= math.floor(math.log10(max))
		log_range			= max_oom - min_oom 
		if log_range	> 1:
			range	= log_range
			min		= min_oom
			max		= max_oom

	# detect order_of_magnitude
	try:
		min_oom					= np.min([math.floor(math.log10(np.abs(min))), math.floor(math.log10(np.abs(max)))])
	except(ValueError):
		if min == 0:
			min_oom					= np.min([0, math.floor(math.log10(np.abs(max)))])
		else:
			min_oom					= np.min([math.floor(math.log10(np.abs(min))), 0])

	detected_oom			= math.floor(math.log10(range))
	if detected_oom == min_oom + 1:
		detected_oom = min_oom

	detected_oom = float(detected_oom)

	# determine spacing
	spacing						= math.floor(range/((number_of_ticks-1) * 10**detected_oom)) * 10**detected_oom
	if spacing == 0:
		spacing					= math.ceil(range/((number_of_ticks-1) * 10**detected_oom)) * 10**detected_oom
	
	# determine ticks
	tick_min					= min - min%spacing
	tick_max					= max - max%spacing
	if ("log" in scale) and (log_range > 1):
		ticks							= np.logspace(tick_min, tick_max, math.ceil((tick_max - tick_min)/spacing)+1)
	else:
		ticks							= np.linspace(tick_min, tick_max, math.ceil((tick_max - tick_min)/spacing)+1)
	if ("log" in scale) and (ticks[0] < min/10):
		ticks						= ticks[1:]		

	# print("\n\n**********\n\n", min, "\n", max, "\n", ticks, "\n\n**********\n\n")

	# set ticklabels
	if tickformat == "power_of_10":
		ticklabels	= []
		for tick in ticks:
			ticklabels.append("${}$".format(ticker.ScalarFormatter(useOffset=False, useMathText=True)._formatSciNotation("%1.10e" %(tick))))
	else:
		ticklabels	= [tickformat %(tick) for tick in ticks]
	
	return ticks, ticklabels


############################################

'''
##########################################
##                                                                                         ##
##                            General Results Plotting                           ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to plot_python_analysis_results     ##########

def plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color, z=[], ):

	
    if "sequence" in sample_plot_details.action:
        
        x = sample_gudrun_results.parameter
        xlabel = "Parameter"                    # To be changed by user
        x_width = 0.5                           # width for color bars in color scale plots
        filename += "Parameter"                 # To be changed by user
    
    
        
    else:
        
        x = sample_gudrun_results.sample
        xlabel = "Sample"
        x_width	= 5.				
        filename += "Sample"

	# Margins & positions
    top				= 0.92
    bottom			= 0.14
    left				= 0.14
    if plot == "colormap":
    	right			= 0.88
    	hspace		= 0.02
    	cbarwidth	= 0.03
    else:
    	right			= 0.998
    	hspace		= 0.
    	cbarwidth	= 0.

	# Initialise Plot:
    fig = plt.figure(figsize=(9,6))
    fig.patch.set_facecolor("w")					# colour of outer box
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)	# margins
    mpl.rcParams["lines.linewidth"]	= 2			# linewidth of plots
    plt.rcParams["mathtext.default"] = "regular"		# mathtext same font as regular text
	#~ rc('axes', linewidth=2)						# linewidth of axes
    mpl.rcParams["axes.linewidth"]	= 2

	# Subplot positions
	# [x position, y position, rel width, rel heigth]
    ax		= plt.axes([left, bottom, right-(cbarwidth+hspace+left), top-bottom])	
	# colorbar
    if plot == "colormap":
    		cbar			= plt.axes([right-cbarwidth, bottom, cbarwidth, top-bottom])	

	# plot data
    if plot=="fill_between":
    	ax.fill_between(x, y[0], y[1], color='k', zorder=0)
    elif plot=="errorbar":
    	ax.errorbar(x, y[0],  yerr=y[1], capthick=2, capsize=5, color='k', zorder=0)
    elif plot=="scatter_errorbar":

        ax.errorbar(x, y[0],  yerr=y[1], capthick=2, capsize=5, color='k', zorder=0)
        ax.plot(x, y[0], color='k', marker="o", markersize=8, zorder=0)

    elif plot == "colormap":
       # pcolor requires x and y to be boundaries for z value field, rather than data point at which to plot z-values
    	# therefore, create new arrays X and Y that will contain boundaries rather than data points (plus increased dimension)
		# example: x = [1, 2, 3]		=> 	X = [[0.5, 1.5, 2.5, 3.5], [0.5, 1.5, 2.5, 3.5]] 	=> 	z(1) will be plotted across x = 0.5 to x = 1.5

		# handling data in subsets will adjust color scale to match range of any given subset
		# therefore, append max and min values of total data set to each subset (using dummy x and y values at plotting range limits)

		# plot ranges
        z_min = 0.
        z_max = 100.
    				
		# set z values
        z_expanded = []
        
        for i in range(len(z)):
            z_expanded.append([z_min])
            
            for j in range(len(z[i])):
                z_expanded[i].append(z[i][j])
                
        z_expanded[i].append(z_max)
        z_expanded[i].append(z_max)
            
      
		# set x borders
        for i in range(len(x)):
            if i > 0 :
                x_left	= (x[i-1]+x[i])/2.
            else:
                try:
                    x_left	= np.min([x[i]-(x[i]-x[i+1])/2., x[i]-x_width])
                except(IndexError):
                    x_left	= x[i]-x_width
                
            try:
                x_right =(x[i]+x[i+1])/2.
            except(IndexError):
                x_right = np.max([x[i]+(x[i-1]-x[i])/2., x[i]+x_width])
                

			# set y borders
            y_expanded	= np.zeros(len(y[i])+3)
            for j in range(len(y[i])):
                if j == 0:
                    if len(y[i]) > 1:
                        y_expanded[j]		= y[i][j]-(y[i][j]-y[i][j+1])/2.
                        y_expanded[j+1]		= y[i][j]+(y[i][j]-y[i][j+1])/2
                        
                    else:
                        y_expanded[j]		= y[i][j]-0.01
                        y_expanded[j+1]	= y[i][j]+0.01
                        
                else:
                    y_expanded[j+1]		= (y[i][j-1]+y[i][j])/2.
                    
            if len(y[i]) > 1:
                y_expanded[-2]	= y[i][-1]+ (y[i][-1] - y[i][-2])/2.
                y_expanded[-1]	= y[i][-1]+ (y[i][-1] - y[i][-2])/2.
            y_expanded		= np.array(y_expanded)
            
            # set arrays
            
            X = np.transpose([np.zeros(len(y[i])+3) + x_left, np.zeros(len(y[i])+3) + x_right])
            Y = np.transpose([y_expanded, y_expanded])
            Z = np.transpose([z_expanded[i]])
            
            # Plot
        
            ax.pcolormesh(X, Y, Z, cmap=color[0])
            
        # colorbar
            
        cb = mpl.colorbar.ColorbarBase(
					cbar											, 
					cmap			= color[0]				,
					norm				= color[1]				,
					boundaries	= color[2]				,
					ticks				= color[3]				,
					spacing			= 'proportional'	,
					orientation		= 'vertical'
					)
        
        # ticks
        cbar.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
        #label
        cbar.text(1.0,bottom+(top-bottom)/2., color[4], transform=fig.transFigure, horizontalalignment='right', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
        
    # labels
    
    # xlabel
    ax.text(left+(right-left)/2.,0., xlabel, transform=fig.transFigure, horizontalalignment='center', verticalalignment='bottom', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
    # ylabel
    ax.text(0.0, bottom+(top-bottom)/2., ylabel, transform=fig.transFigure, horizontalalignment='left', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
    
    # axes
    if "Parameter" in xlabel:          
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmin, xmax)
        
    
    ax.set_ylim(yscale[0], yscale[1])
    
    try:
        if "log" in yscale[2]:
            ax.set_yscale('log')
        ax.yaxis.set_ticks(yscale[3])
        ax.yaxis.set_ticklabels(yscale[4])
    except(IndexError):
        pass
    
    # tick appearances
    ax.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
    ax.tick_params(which="minor", labelsize=sample_plot_details.minor_tick_fontsize, width=sample_plot_details.minor_tick_width, length=sample_plot_details.minor_tick_length, direction="out", labelleft=False)
    
    # save plot
    plt.savefig(filename + ".pdf")
    plt.savefig(filename + ".png")
    plt.close()
    
    return


############################################
