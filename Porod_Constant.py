import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib.colors import *
from matplotlib import animation
import time
from copy import *
from analysis_input_file_handling import *
from plot_python_analysis_results import *
from plot_settings import *


'''
##########################################
##                                                                                         ##
##                      DATA & PARAMETER HANDLING                      ##
##                                                                                         ##
##########################################
'''
	

###########     FUNCTION to set_porod_constant_fit_ranges     ##########

def set_porod_constant_fit_ranges(sample_analysis_directories, sample_plot_details, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results):
    
    if "sequence" in sample_plot_details.action:
        
        os.chdir(sample_analysis_directories.Python_directory)
        parameter, Porod_range_left, Porod_range_right	= np.loadtxt("GP_analysis_Porod_range_%s_%s.txt" %(sample_plot_details.sample, sample_plot_details.action), unpack=True)
            
    else:
        
        os.chdir(sample_analysis_directories.Python_directory)
        sample, Porod_range_left, Porod_range_right	= np.loadtxt("GP_analysis_Porod_range_%s_%s.txt" %(sample_plot_details.sample, sample_plot_details.action), unpack=True)
        
    
    for i in range(sample_porod_constant_results.number_of_files):
        
        file = sample_gudrun_results.file_start[i]
        
        if file in sample_porod_constant_input.files_manual_fit_range:
            pass
        
        else:
            try:
                sample_porod_constant_input.q_fit_range_min[file]		= Porod_range_left[i] * (1. - sample_porod_constant_input.expand_fit_range)
                sample_porod_constant_input.q_fit_range_max[file]	= Porod_range_right[i] * (1. + sample_porod_constant_input.expand_fit_range)
            except(IndexError):
                sample_porod_constant_input.q_fit_range_min[file]		= Porod_range_left * (1. - sample_porod_constant_input.expand_fit_range)
                sample_porod_constant_input.q_fit_range_max[file]	= Porod_range_right * (1. + sample_porod_constant_input.expand_fit_range)
                
    return
            

	

############################################


'''
##########################################
##                                                                                         ##
##                                     PLOTTING                                      ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to plot_individual_porod_constant_fits     ##########

def plot_individual_porod_constant_fits(i, sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results):
    
    
    if "sequence" in sample_plot_details.action:
        
            
        file_start	= sample_gudrun_results.file_start[i]
        parameter = sample_gudrun_results.parameter[i]
    
        q_min		= sample_porod_constant_input.q_fit_range_min[file_start]
        q_max		= sample_porod_constant_input.q_fit_range_max[file_start]
    
        q = sample_porod_constant_results.q_data[i]
        y = sample_porod_constant_results.y_data[i]
        yerr = sample_porod_constant_results.y_data_err[i]
        q_fit_start	= sample_porod_constant_results.q_fit_start[i]
        q_fit_end	= sample_porod_constant_results.q_fit_end[i]
        K_line		= sample_porod_constant_results.K_line[i]
    
        K = sample_porod_constant_results.K[i]
        K_err = sample_porod_constant_results.K_err[i]
        SSA = sample_porod_constant_results.SSA[i]
        SSA_err	= sample_porod_constant_results.SSA_err[i]
    
    
        
    else:
    
        file_start	= sample_gudrun_results.file_start[i]
        sample = sample_gudrun_results.sample[i]
    
        q_min		= sample_porod_constant_input.q_fit_range_min[file_start]
        q_max		= sample_porod_constant_input.q_fit_range_max[file_start]
    
        q = sample_porod_constant_results.q_data[i]
        y = sample_porod_constant_results.y_data[i]
        yerr = sample_porod_constant_results.y_data_err[i]
        q_fit_start	= sample_porod_constant_results.q_fit_start[i]
        q_fit_end	= sample_porod_constant_results.q_fit_end[i]
        K_line		= sample_porod_constant_results.K_line[i]
    
        K = sample_porod_constant_results.K[i]
        K_err = sample_porod_constant_results.K_err[i]
        SSA = sample_porod_constant_results.SSA[i]
        SSA_err	= sample_porod_constant_results.SSA_err[i]
    
    # labels
    
    xlabel = "Q ($\AA^{-1}$)"
    ylabel = "I(Q)$\cdot$Q$^4$ (m$^2$/cm$^3$/$\AA^4$)"
    
    
    if "sequence" in sample_plot_details.action:
            
        plot_title = "Scan %i, %.1f " %(file_start, parameter)
                
    else:
        
        plot_title = "Scan %i, %i" %(file_start, sample)
 

        
    # margins and positions
    
    top = 0.92
    bottom = 0.14
    left = 0.12
    right = 0.83
    
    # Initialise Plot:
    
    fig = plt.figure(figsize=(9,6))
    fig.patch.set_facecolor("w")					# colour of outer box
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)	# margins
    mpl.rcParams["lines.linewidth"]	= 2			# linewidth of plots
    plt.rcParams["mathtext.default"] = "regular"		# mathtext same font as regular text
    rc('axes', linewidth=2)						# linewidth of axes
    ax = plt.subplot(1,1,1)
    
    # plot fit
    
    ax.plot(q, K_line, 'k-')
    ax.plot(q, K_line+K_err,'k--')
    ax.plot(q, K_line-K_err, 'k--')
    result_parameters_1	= 		"SSA\n\n\nK\n"
    result_parameters_2	= 		"       = %.2f\n       $\pm$ %.2f\n\n       = %.1e\n       $\pm$ %.1e" %(SSA, SSA_err, K, K_err)
    
    # plot data
    
    ax.plot(q, y, color='yellow', linewidth=2)
    ax.fill_between(q, y-yerr, y+yerr, color='yellow', alpha=0.5, zorder=2)
    
    # axes
    
    xmin = 0
    xmax = np.max(q) + (np.max(q) - np.min(q))*0.05
    ymin = 0
    ymax = np.max(y) + (np.max(y) - np.min(y))*0.1
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    
    #labels
    
    # title
    ax.text(0.5, 1.0, plot_title, transform=fig.transFigure, horizontalalignment='center', verticalalignment='top', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
    
    ax.text(left+(right-left)/2.,0., xlabel, transform=fig.transFigure, horizontalalignment='center', verticalalignment='bottom', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
    
    ax.text(0.0,bottom+(top-bottom)/2., ylabel, transform=fig.transFigure, horizontalalignment='left', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
    
    # results parameters
    ax.text(right*1.01, bottom+(top-bottom)*0.97, result_parameters_1, transform=fig.transFigure, horizontalalignment="left", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
    ax.text(right*1.01, bottom+(top-bottom)*0.97, result_parameters_2, transform=fig.transFigure, horizontalalignment="left", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
    # points of interest
    ax.text((q_fit_start + q_fit_end)/2., ymin + (ymax-ymin)*0.275, "fit\nrange", horizontalalignment="center", verticalalignment="center", fontsize=sample_plot_details.annotation_fontsize)
    ax.text((q_min + q_max)/2., ymin + (ymax-ymin)*0.1, "range\nconsidered\nfor fit", horizontalalignment="center", verticalalignment="center", fontsize=sample_plot_details.annotation_fontsize)
    
    # shaded area to indicate fit range
    ax.fill_between([q_fit_start, q_fit_end], [ymin + (ymax-ymin)*0.2, ymin + (ymax-ymin)*0.2], [ymax, ymax], color="0.8", zorder=1)
    ax.fill_between([q_min, q_max], [ymin, ymin], [ymax, ymax], color="0.9", zorder=0)
    
    #tick appearances
    ax.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
    ax.tick_params(which="minor", labelsize=sample_plot_details.minor_tick_fontsize, width=sample_plot_details.minor_tick_width, length=sample_plot_details.minor_tick_length, direction="out", labelleft=False)
    
    #save plot
    
    plt.savefig("Porod_constant_NIMROD000%i.png" %(file_start))
    plt.close()
    
    return
    


###########     FUNCTION to plot_porod_range_development_in_Porod_Constant_fits     ##########

def plot_porod_range_development_in_Porod_Constant_fits(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results):

	# set parameters for plot
	filename	= 	"%s_%s__Porod_range_vs_" %(sample_plot_details.sample, sample_plot_details.action)
	ylabel		= "Porod Fit Range Q ($\AA^{-1}$)"
	# color			[color, alpha]
	color			= ["k", 0.5]
	# axes			[ymin, ymax, logscale?, tickpositions?, ticklables?]
	yscale		= [sample_porod_constant_input.q_plot_min, sample_porod_constant_input.q_plot_max, "log", [0.02, 0.04, 0.07, 0.1, 0.2], ["0.02", "0.04", "0.07", "0.10", "0.20"]]

	# plot settings and y-data
	if len(sample_porod_constant_results.q_fit_start)>1:
		plot	= "fill_between"
		# [lower range end, upper range end]
		y		= [sample_porod_constant_results.q_fit_start, sample_porod_constant_results.q_fit_end]
	else:
		plot	= "errorbar"
		# [y values, y errors]
		y 		= [	(sample_porod_constant_results.q_fit_start + sample_porod_constant_results.q_fit_end) / 2.	,
						(sample_porod_constant_results.q_fit_end - sample_porod_constant_results.q_fit_start) / 2.	]
					
	plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
	
	return


###########     FUNCTION to plot SSA     ##########

def plot_Porod_constant_SSA(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results):

	# set parameters for plot
	filename	= 	"%s_%s__SSA_vs_" %(sample_plot_details.sample, sample_plot_details.action)
	ylabel		= "Specific Surface Area (m$^2$/cm$^3$)"
	# color			[color, alpha]
	color			= ["k", 0.5]
	# axes			[ymin, ymax, logscale?, tickpositions?, ticklables?]
	yscale		= [0, 150.]     # To be changed by user

	SSA			= sample_porod_constant_results.SSA
	SSA_err	= sample_porod_constant_results.SSA_err

	# plot settings and y-data
	plot	= "scatter_errorbar"
	# [[y values, y errors], [ , ], ...]
	y 		= [SSA, SSA_err]	
					
	plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
	
	# logscale
	yscale = [1., 1000]     # To be changed by user                         
	yscale.append("log")
	filename	= filename.split("_vs_")[0] + "_logscale_vs_"
	plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
	
	return


###########     FUNCTION to plot all porod constants     ##########

def plot_all_Porod_constants(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results):

    # set parameters for plot
    filename    =   "%s_%s__Porod_constant" %(sample_plot_details.sample, sample_plot_details.action)
    xlabel      = "Q ($\AA^{-1}$)"
    ylabel      = "I(Q)$\cdot$Q$^4$ (m$^2$/cm$^3$/$\AA^4$)"

    # Margins & positions
    top                 = 0.92
    bottom              = 0.14
    left                    = 0.15
    right                   = 0.88
    hspace              = 0.02
    cbarwidth           = 0.03


    if "sequence" in sample_plot_details.action:

        # exception handling required for cases with only one data point
        try:
            n   = len(sample_gudrun_results.file_start)
        except (TypeError):
            n    = 1


        spectrum = plt.get_cmap('gnuplot')
        colors = spectrum(np.linspace(0, 1, sample_porod_constant_results.number_of_files))

    
        # Initialise Plot:

        fig= plt.figure(figsize=(9,6))
        fig.patch.set_facecolor("w")                                                         # colour of outer box
        plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)  # margins
        mpl.rcParams["lines.linewidth"]  = 2                                         # linewidth of plots
        plt.rcParams["mathtext.default"] = "regular"                                 # mathtext same font as regular text
        rc('axes', linewidth=2)

        # Subplot positions
        # [x position, y position, rel width, rel heigth]
        ax_DCS       = plt.axes([left, bottom, right-(cbarwidth+hspace+left), top-bottom])   
        # colorbar
        cbar         = plt.axes([right-cbarwidth, bottom, cbarwidth, top-bottom])

        for i in range(n):

            if sample_gudrun_results.file_start[i] in sample_plot_details.exceptions:
                print ("NIMROD000%i.mint01 not plotted" %(sample_gudrun_results.file_start[i]))

            else:
                ax_DCS.plot(sample_porod_constant_results.q_data[i], sample_porod_constant_results.y_data[i], color=colors[i])
                ax_DCS.plot(sample_porod_constant_results.q_data[i], sample_porod_constant_results.K_line[i], '--', color=colors[i])

        # colorbar
        cb = mpl.colorbar.ColorbarBase(
                cbar,
                cmap = plt.cm.gnuplot,
                norm = sample_plot_details.color_norm,
                boundaries = sample_plot_details.color_bounds,
                ticks = sample_plot_details.color_plot_ticks,
                spacing = 'proportional',
                orientation = 'vertical'
                )

        # tick appearances

        ax_DCS.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
        ax_DCS.tick_params(which="minor", labelsize=sample_plot_details.minor_tick_fontsize, width=sample_plot_details.minor_tick_width, length=sample_plot_details.minor_tick_length, direction="out", labelleft=False)
        #colorbar
        cbar.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")

        #labels

        ax_DCS.text(0.5, 1.0, sample_plot_details.plot_title, transform=fig.transFigure, horizontalalignment='center', verticalalignment='top', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
        # xlabel
        ax_DCS.text(left+(right-left)/2., 0., xlabel, transform=fig.transFigure, horizontalalignment='center', verticalalignment='bottom', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
        # ylabel
        ax_DCS.text(0.0, bottom+(top-bottom)/2., ylabel, transform=fig.transFigure, horizontalalignment='left', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
        # colorbar
        cbar.text(1.0, bottom+(top-bottom)/2., sample_plot_details.color_bar_label, transform=fig.transFigure, horizontalalignment='right', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)

        # axes

        xmin, xmax   = ax_DCS.get_xlim()
        ymin, ymax   = ax_DCS.get_ylim()
        ax_DCS.set_xlim(0., xmax)
        ax_DCS.set_ylim(0., ymax)

        # save plot

        plt.savefig(filename + ".pdf")
        plt.savefig(filename + ".png")

        plt.close()

        return


    else:

        # exception handling required for cases with only one data point
        try:
            n   = len(sample_gudrun_results.file_start)
        except (TypeError):
            n    = 1


        spectrum = plt.get_cmap('gnuplot')
        colors = spectrum(np.linspace(0, 1, sample_porod_constant_results.number_of_files))

        # Initialise Plot:

        fig= plt.figure(figsize=(9,6))
        fig.patch.set_facecolor("w")                                                         # colour of outer box
        plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)  # margins
        mpl.rcParams["lines.linewidth"]  = 2                                         # linewidth of plots
        plt.rcParams["mathtext.default"] = "regular"                                 # mathtext same font as regular text
        rc('axes', linewidth=2)



        # Subplot positions
        # [x position, y position, rel width, rel heigth]
        ax_DCS       = plt.axes([left, bottom, right-(cbarwidth+hspace+left), top-bottom])   
        
        for i in range(n):

            if sample_gudrun_results.file_start[i] in sample_plot_details.exceptions:
                print ("NIMROD000%i.mint01 not plotted" %(sample_gudrun_results.file_start[i]))

            else:

                ax_DCS.plot(sample_porod_constant_results.q_data[i], sample_porod_constant_results.y_data[i], label=sample_gudrun_results.sample[i], c=colors[i])
                ax_DCS.plot(sample_porod_constant_results.q_data[i], sample_porod_constant_results.K_line[i], '--', c=colors[i])

    
        # tick appearances

        ax_DCS.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
        ax_DCS.tick_params(which="minor", labelsize=sample_plot_details.minor_tick_fontsize, width=sample_plot_details.minor_tick_width, length=sample_plot_details.minor_tick_length, direction="out", labelleft=False)
        
        #labels

        ax_DCS.text(0.5, 1.0, sample_plot_details.plot_title, transform=fig.transFigure, horizontalalignment='center', verticalalignment='top', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
        # xlabel
        ax_DCS.text(left+(right-left)/2., 0., xlabel, transform=fig.transFigure, horizontalalignment='center', verticalalignment='bottom', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
        # ylabel
        ax_DCS.text(0.0, bottom+(top-bottom)/2., ylabel, transform=fig.transFigure, horizontalalignment='left', verticalalignment='center', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
        
        # axes

        xmin, xmax   = ax_DCS.get_xlim()
        ymin, ymax   = ax_DCS.get_ylim()
        ax_DCS.set_xlim(0., xmax)
        ax_DCS.set_ylim(0., ymax)

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, frameon=False)

        # save plot

        plt.savefig(filename + ".pdf")
        plt.savefig(filename + ".png")

        plt.close()

        return
             

         



	  

    






############################################


'''
##########################################
##                                                                                         ##
##                                  RUNNING FITS                                  ##
##                                                                                         ##
##########################################
'''

###########     FUNCTION to run_porod_constant_fit     ##########

def run_porod_constant_fit(i, sample_material_properties, sample_plot_details, sample_porod_constant_temporary_fit_parameters, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results):

	print("Porod Constant fit")

	file_start	= sample_gudrun_results.file_start[i]

	#Convert DCS (barn/sr/atom) to I(q) (m^2/cm^3)		 (barn & angstrom conversion: 1e-28*1e24=1e-4 )
	I 			= (sample_porod_constant_temporary_fit_parameters.y*sample_material_properties.atomic_number_density*1e-4)	
	Ierr		= (sample_porod_constant_temporary_fit_parameters.yerr*sample_material_properties.atomic_number_density*1e-4)

	# set fit parameters
	q				= sample_porod_constant_temporary_fit_parameters.q
	q_min		= sample_porod_constant_input.q_fit_range_min[file_start]
	q_max		= sample_porod_constant_input.q_fit_range_max[file_start]
	fit_points	= sample_porod_constant_input.number_of_points_to_fit_over
	#Convert I(q) to I(q)*q^4 (m^2/cm^3*angstrom^4) for fit:
	y				= I*q**4
	yerr			= Ierr*q**4

	start_index	= np.where(q>=q_min)[0][0]
	stop_index	= np.where(q<=q_max)[0][-1] + 1
	# catch case where Porod range is smaller than given number of fit points:
	if stop_index - start_index < 2*fit_points:
		if stop_index - 2*fit_points - 1 >= 0:
			start_index		= stop_index - 2*fit_points - 1
			sample_porod_constant_input.q_fit_range_min[file_start]	= q[start_index]
		else:
			stop_index		= start_index + 2*fit_points + 1
			sample_porod_constant_input.q_fit_range_max[file_start]	= q[stop_index]

	# Find best fit of constant curve to data
	for l in range(start_index, stop_index - fit_points):
		p0			= 1e-8
		x_fit 		= q[l:l+fit_points]
		y_fit 		= y[l:l+fit_points]
		y_fit_err 	= yerr[l:l+fit_points]
		weight		= y_fit_err/y_fit
		
		result, covariant = curve_fit(lambda x, p0: p0, x_fit, y_fit, p0=p0, sigma=weight, absolute_sigma=True, maxfev=10000)

		y_err	= np.sqrt((y_fit_err**2).sum())/len(y_fit_err)
		if (len(y_fit) > 1) and covariant is not None:
			error	= np.sqrt(((y_fit - result)**2).sum()/(len(y_fit)-1))
		else:
			error = 0
		try:
			result	= result[0]
		except(IndexError):
			pass

		sample_porod_constant_temporary_fit_parameters.fit_attempts	+= 1
		sample_porod_constant_temporary_fit_parameters.q_fit_ranges.append([x_fit[0], x_fit[-1]] )
		sample_porod_constant_temporary_fit_parameters.K.append(result)
		sample_porod_constant_temporary_fit_parameters.K_err.append(error)
		sample_porod_constant_temporary_fit_parameters.y_err.append(y_err)

		#~ print("\t%.2e +/- %.2e" %(result, error))
		# print("\n")
		# print("result\t\t%.2e" %result)
		# print("error\t\t%.2e" %error)
		# print("best fit error\t%.2e" %sample_porod_constant_results.K_err[i])

		if (error < sample_porod_constant_results.K_err[i]) and (result > 5e-10):
			sample_porod_constant_results.K[i]				= result
			sample_porod_constant_results.K_err[i]			= np.max([error, y_err])
			sample_porod_constant_results.q_fit_start[i]	= x_fit[0]
			sample_porod_constant_results.q_fit_end[i]		= x_fit[-1]

	print("Fit result: K = %.2e +/- %.2e" %(sample_porod_constant_results.K[i], sample_porod_constant_results.K_err[i]))

	# write data to lists
	sample_porod_constant_results.q_data.append(np.array(q))
	sample_porod_constant_results.y_data.append(np.array(y))
	sample_porod_constant_results.y_data_err.append(np.array(yerr))
	sample_porod_constant_results.K_line.append(np.array(np.zeros(len(y))+ sample_porod_constant_results.K[i]))
	sample_porod_constant_results.fit_deviations.append(np.array(np.abs(y-sample_porod_constant_results.K[i])/sample_porod_constant_results.K[i]))
	# Compute SSA (m^2/cm^3):
	sample_porod_constant_results.SSA[i]				= sample_porod_constant_results.K[i]/(2.0*np.pi*sample_material_properties.scattering_length_density_difference**2.0)
	sample_porod_constant_results.SSA_err[i]			= sample_porod_constant_results.K_err[i]/(2.0*np.pi*sample_material_properties.scattering_length_density_difference**2.0)

	return
	
	
###########     FUNCTION to call_porod_constant_analysis_loop     ##########

def call_porod_constant_analysis_loop(i, sample_material_properties, sample_analysis_directories, sample_porod_constant_temporary_fit_parameters, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results, sample_plot_details):

    if "sequence" in sample_plot_details.action:

        file_start= sample_gudrun_results.file_start[i]
        parameter               = sample_gudrun_results.parameter[i]
        q_fit_min   = sample_porod_constant_input.q_fit_range_min[file_start]
        q_fit_max   = sample_porod_constant_input.q_fit_range_max[file_start]

    else:

        file_start= sample_gudrun_results.file_start[i]
        sample               = sample_gudrun_results.sample[i]
        q_fit_min   = sample_porod_constant_input.q_fit_range_min[file_start]
        q_fit_max   = sample_porod_constant_input.q_fit_range_max[file_start]

    # read neutron data from file
    python_fits_read_neutron_data_from_file(i, [q_fit_min, q_fit_max], sample_analysis_directories, sample_porod_constant_input, sample_gudrun_results, sample_plot_details, sample_porod_constant_temporary_fit_parameters, sample_porod_constant_results)

    # Fit data:

    if file_start in sample_plot_details.exceptions:
        sample_porod_constant_results.q_data.append(np.array([0]))
        sample_porod_constant_results.y_data.append(np.array([0]))
        sample_porod_constant_results.y_data_err.append(np.array([0]))
        sample_porod_constant_results.fit_deviations.append(np.array([0]))
        sample_porod_constant_results.K_line.append(np.array([0]))

    else:
        sample_porod_constant_results.K_err[i]                              = 1e100
        popt = run_porod_constant_fit(i, sample_material_properties, sample_plot_details, sample_porod_constant_temporary_fit_parameters, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results)

        # Plot fit:

        plot_individual_porod_constant_fits(i, sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results)
        
    return


############################################

'''
##########################################
##                                                                                         ##
##                                SAVING RESULTS                                ##
##                                                                                         ##
##########################################
'''


#############     FUNCTION to save results to file     ############

def save_porod_constant_results_to_file(input, sample_plot_details, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results):
    
    if "sequence" in sample_plot_details.action:
            
        Save_Header = ""
        Save_Data = []
        Save_Format = []
            
        # sort results into variable
            
        Save_Header += "parameter, "
        Save_Data.append(sample_gudrun_results.parameter)
        Save_Format.append('%1.4e')
        Save_Header += "SSA, "
        Save_Data.append(sample_porod_constant_results.SSA)
        Save_Format.append('%1.4e')
        Save_Header		+= "SSA_err, "
        Save_Data.append(sample_porod_constant_results.SSA_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "K, "
        Save_Data.append(sample_porod_constant_results.K)
        Save_Format.append('%1.4e')
        Save_Header		+= "K_err, "
        Save_Data.append(sample_porod_constant_results.K_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "q_fit_start, "
        Save_Data.append(sample_porod_constant_results.q_fit_start)
        Save_Format.append('%1.4e')
        Save_Header		+= "q_fit_end, "
        Save_Data.append(sample_porod_constant_results.q_fit_end)
        Save_Format.append('%1.4e')
            
        # save results to file
        np.savetxt("SSA_+_Porod_analysis_%s_%s.txt" %(input["sample"], input["action"]), np.transpose(Save_Data), fmt=Save_Format, header=Save_Header)
            
    else:
        
        Save_Header = ""
        Save_Data = []
        Save_Format = []
            
        # sort results into variable
            
        Save_Header += "sample, "
        Save_Data.append(sample_gudrun_results.sample)
        Save_Format.append('%1.4e')
        Save_Header += "SSA, "
        Save_Data.append(sample_porod_constant_results.SSA)
        Save_Format.append('%1.4e')
        Save_Header		+= "SSA_err, "
        Save_Data.append(sample_porod_constant_results.SSA_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "K, "
        Save_Data.append(sample_porod_constant_results.K)
        Save_Format.append('%1.4e')
        Save_Header		+= "K_err, "
        Save_Data.append(sample_porod_constant_results.K_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "q_fit_start, "
        Save_Data.append(sample_porod_constant_results.q_fit_start)
        Save_Format.append('%1.4e')
        Save_Header		+= "q_fit_end, "
        Save_Data.append(sample_porod_constant_results.q_fit_end)
        Save_Format.append('%1.4e')
            
        # save results to file
        np.savetxt("SSA_+_Porod_analysis_%s_%s.txt" %(input["sample"], input["action"]), np.transpose(Save_Data), fmt=Save_Format, header=Save_Header)
            
            

############################################

'''
##########################################
##                                                                                         ##
##                                 MAIN FUNCTION                                ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to run_porod_constant_analysis     ##########

def run_porod_constant_analysis(input, sample_material_properties, sample_analysis_directories, sample_plot_details, sample_gudrun_results, sample_porod_constant_input):
	print("\n\n-------------------\n\nRunning Porod Constant analysis.\n\n")

	# Initialise arrays for fit results
	try:
		n	= len(sample_gudrun_results.file_start)
	except(TypeError):
		n	= 1
	sample_porod_constant_results = porod_constant_results(n)
		
	# Set q fit range (based on double Guinier-Porod fit):
	set_porod_constant_fit_ranges(sample_analysis_directories, sample_plot_details, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results)	

	for i in range(sample_porod_constant_results.number_of_files):
		# Initialise temporary storage for data and start parameters
		sample_porod_constant_temporary_fit_parameters = porod_constant_temporary_fit_parameters()
		# Call fit
		call_porod_constant_analysis_loop(i, sample_material_properties, sample_analysis_directories, sample_porod_constant_temporary_fit_parameters, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results, sample_plot_details)

	print("----------------------------------\n----------------------------------\n\nFits complete\n")

	# Save results to file
	save_porod_constant_results_to_file(input, sample_plot_details, sample_porod_constant_input, sample_gudrun_results, sample_porod_constant_results)
	
	# Plot results
	plot_Porod_constant_SSA(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results)
	plot_porod_range_development_in_Porod_Constant_fits(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results)
	plot_all_Porod_constants(sample_plot_details, sample_porod_constant_input, sample_porod_constant_results, sample_gudrun_results)

