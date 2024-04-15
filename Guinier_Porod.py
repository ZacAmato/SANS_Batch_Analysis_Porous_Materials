import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
import math
from copy import *
from analysis_input_file_handling import *
from plot_python_analysis_results import *
from plot_settings import *

'''
##########################################
##                                                                                         ##
##                                 FIT FUNCTIONS                                 ##
##                                                                                         ##
##########################################
'''



#############     FUNCTION guinier_porod     ############

def guinier_porod_function(x, G, s, Rg, d, q_min, q_max):
	if ((d-s)*(3.0-s) > 0) and (Rg > 2./q_max) and (s > 0) and (d > 0):
		Q1	= (1.0/Rg)*np.sqrt(((d-s)*(3.0-s))/2.0)
		A	= G * Q1**(d-s) * np.exp(-Q1**2.*Rg**2./(3.-s))
		if (q_min<Q1) and (Q1 < q_max) and (A>0) and (G>0):
			GP	= np.piecewise(x, [x <= Q1, x > Q1], [lambda x: (G/x**s)*np.exp((-x**2.0*Rg**2.0)/(3.0-s)), lambda x: A/x**d])
			result = np.log(GP)
			# print("fitting GP")
		else:
			result = np.zeros(np.size(x)) + 1e10
			print("Q1 out of range")
	else:
		result = np.zeros(np.size(x)) + 1e10
		if d <= 0:
			print("\nWARNING\nd > 0:\t%s" %(d>0), d)
		if (d-s)*(3.-s) <= 0:
			print("\nWARNING\n(d-s)*(3.-s) > 0:\t%s" %((d-s)*(3.-s) > 0), d, ", ", s)
		if Rg <= 2./q_max:
			print("\nWARNING\nRg > 2/q_max:\t%s" %(Rg>2./q_max) , Rg, ", ", 2./q_max)
		if s <= 0:
			print("\nWARNING\ns > 0:\t%s" %(s>0), s)
	return result



############################################

'''
##########################################
##                                                                                         ##
##                           PARAMETER HANDLING                            ##
##                                                                                         ##
##########################################
'''

def assign_manual_guinier_porod_fit_parameters(i, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results):
	file_start						= sample_gudrun_results.file_start[i]
	uncertainty						= sample_guinier_porod_input.manual_fit_uncertainty[file_start]/100.
	#~ sample_guinier_porod_results.fit_tag[i]			= 0
	
	sample_guinier_porod_results.d[i]					= sample_guinier_porod_input.manual_fit_d[file_start]
	sample_guinier_porod_results.G[i]				= sample_guinier_porod_input.manual_fit_G_low_Q[file_start]
	sample_guinier_porod_results.d_err[i]			= uncertainty * sample_guinier_porod_input.manual_fit_d[file_start]
	sample_guinier_porod_results.G_err[i]		= uncertainty * sample_guinier_porod_input.manual_fit_G_low_Q[file_start]
		
	sample_guinier_porod_results.s[i]			= sample_guinier_porod_input.manual_fit_s_low_Q[file_start]
	sample_guinier_porod_results.Rg[i]		= sample_guinier_porod_input.manual_fit_Rg_low_Q[file_start]
	sample_guinier_porod_results.s_err[i]	= uncertainty * sample_guinier_porod_input.manual_fit_s_low_Q[file_start]
	sample_guinier_porod_results.Rg_err[i]	= uncertainty * sample_guinier_porod_input.manual_fit_Rg_low_Q[file_start]

	sample_guinier_porod_results.porod_range[0][i]	= sample_guinier_porod_input.q_plot_min
	sample_guinier_porod_results.porod_range[1][i]	= sample_guinier_porod_input.q_plot_max
		
	
	
	return 



###########     FUNCTION to assign manual start parameters     ##########

def assign_manual_start_parameters(file_start, sample_guinier_porod_input, sample_guinier_porod_temporary_fit_parameters):
    
        
    sample_guinier_porod_temporary_fit_parameters.p0_porod		= [sample_guinier_porod_input.manual_start_G_low_Q[file_start], sample_guinier_porod_input.manual_start_d[file_start]]
        
    q = 0.02
    s1_start = sample_guinier_porod_input.manual_start_s_low_Q[file_start]
    for i in range(len(sample_guinier_porod_temporary_fit_parameters.q_fit)):
        if sample_guinier_porod_temporary_fit_parameters.q_fit[i] >= q:
            q = sample_guinier_porod_temporary_fit_parameters.q_fit[i]
            I_of_q	= sample_guinier_porod_temporary_fit_parameters.y_fit[i]
    G1_start	= I_of_q * q **s1_start
    sample_guinier_porod_temporary_fit_parameters.p0_normal	= [G1_start, 
				sample_guinier_porod_input.manual_start_s_low_Q[file_start],
				sample_guinier_porod_input.manual_start_Rg_low_Q[file_start],
				sample_guinier_porod_input.manual_start_d[file_start]
				]

        
    return
					


###########     FUNCTION to get_short_names_of_guinier_porod_parameters     ##########

def get_short_names_of_guinier_porod_parameters(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results):
    
    if "sequence" in sample_plot_details.action:
                    
        file_start = sample_gudrun_results.file_start[i]
        parameter = sample_gudrun_results.parameter[i]
            
        G1 = sample_guinier_porod_results.G[i]
        s1 = sample_guinier_porod_results.s[i]
        Rg1 = sample_guinier_porod_results.Rg[i]
        G1_err = sample_guinier_porod_results.G_err[i]
        s1_err = sample_guinier_porod_results.s_err[i]
        Rg1_err = sample_guinier_porod_results.Rg_err[i]
            
              
        d = sample_guinier_porod_results.d[i]
        d_err = sample_guinier_porod_results.d_err[i]
        
    
        q_fit = sample_guinier_porod_temporary_fit_parameters.q_fit
        q = sample_guinier_porod_temporary_fit_parameters.q
        y = sample_guinier_porod_temporary_fit_parameters.y
        yerr = sample_guinier_porod_temporary_fit_parameters.yerr
            
        q_min = sample_guinier_porod_input.q_plot_min
        q_max = sample_guinier_porod_input.q_plot_max
        
            
        return file_start, parameter, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max
        
    else:
        
        file_start = sample_gudrun_results.file_start[i]
        sample = sample_gudrun_results.sample[i]
            
        G1 = sample_guinier_porod_results.G[i]
        s1 = sample_guinier_porod_results.s[i]
        Rg1 = sample_guinier_porod_results.Rg[i]
        G1_err = sample_guinier_porod_results.G_err[i]
        s1_err = sample_guinier_porod_results.s_err[i]
        Rg1_err = sample_guinier_porod_results.Rg_err[i]
            
            
        d = sample_guinier_porod_results.d[i]
        d_err = sample_guinier_porod_results.d_err[i]
        
            
        q_fit = sample_guinier_porod_temporary_fit_parameters.q_fit
        q = sample_guinier_porod_temporary_fit_parameters.q
        y = sample_guinier_porod_temporary_fit_parameters.y
        yerr = sample_guinier_porod_temporary_fit_parameters.yerr
            
        q_min = sample_guinier_porod_input.q_plot_min
        q_max = sample_guinier_porod_input.q_plot_max
            
            
        return file_start, sample, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max
            
            
	

############################################

	
'''
##########################################
##                                                                                         ##
##                                      PLOTTING                                     ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to show guinier-porod fit     ##########

def show_guinier_porod_fit_plot(i, sample_gudrun_results, sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_results):
	
	# get_short_names_of_guinier_porod_parameters
    
    if "sequence" in sample_plot_details.action:
            
        file_start, parameter, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max = get_short_names_of_guinier_porod_parameters(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
            
    else:
        
        file_start, sample, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max = get_short_names_of_guinier_porod_parameters(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
        
    
    # plot
    
    plt.plot(q_fit, np.exp(guinier_porod_function(q_fit, G1, s1, Rg1, d, q_min)), color="k", alpha=0.3, linewidth=2)
    plt.plot(q, y, color='yellow', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("%i") %(file_start)
    plt.show()
    plt.close()
    
    return
            
            


###########     FUNCTION to plot_guinier_porod_fit     ##########

def plot_individual_guinier_porod_fits(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results):

	# get_short_names_of_guinier_porod_parameters
    
    if "sequence" in sample_plot_details.action:
            
        file_start, parameter, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max = get_short_names_of_guinier_porod_parameters(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
            
    else:
        
        file_start, sample, G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err, q_fit, q, y, yerr, q_min, q_max = get_short_names_of_guinier_porod_parameters(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
        
        
    # Labels
    
    xlabel = "Q ($\AA^{-1}$)"
    ylabel = "Differential Cross Section (barns/sr/atom)"
    plot_title = "Scan %i         " %(file_start)
    
    # margins and positions
    
    top = 0.92
    bottom = 0.14
    left = 0.12
    right = 0.83
    
    # plot

    fig, ax = plt.subplots(1, 1, figsize=(9,6))
    fig.subplots_adjust(wspace=0, hspace=0)    # Remove horizontal space between axes
    fig.patch.set_facecolor("w")                    # colour of outer box
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom) # margins
    mpl.rcParams["lines.linewidth"] = 2         # linewidth of plots
    plt.rcParams["mathtext.default"] = "regular"        # mathtext same font as regular text
    rc('axes', linewidth=2)                     # linewidth of axes
        
        
    ax.plot(q_fit, np.exp(guinier_porod_function(q_fit, G1, s1, Rg1, d, q_min, q_max)), color="k", alpha=0.3, linewidth=2)
    results_parameters_3 = "G\n\ns\n\nRg\n\n\nd\n"
    results_parameters_4 = "       = %.2f\n       $\pm$ %.2f\n       = %.2f\n       $\pm$ %.2f\n       = %.2f\n       $\pm$ %.2f\n\n       = %.2f\n       $\pm$ %.2f" %(G1, G1_err, s1, s1_err, Rg1, Rg1_err, d, d_err)
        
    # plot data
    
    ax.plot(q, y, color='yellow', linewidth=2)
    ax.fill_between(q, y-yerr, y+yerr, color='yellow', alpha=0.3, zorder=0)
        
    # labels
        
    ax.text(left+(right-left)/2.,0., xlabel, transform=fig.transFigure, horizontalalignment='center', verticalalignment='bottom', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize)
        
    ax.text(0.0,1., ylabel, transform=fig.transFigure, horizontalalignment='left', verticalalignment='top', multialignment="center", color="k", fontsize=sample_plot_details.label_fontsize, rotation=90)
        
    # manual or forced fit?
        
    if sample_guinier_porod_results.fit_tag[i] == 0:
        ax.text(0.97, 0.9, "manual fit", transform=ax.transAxes, horizontalalignment="right", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
    if sample_guinier_porod_results.forced_fit[i] == 1:
        ax.text(0.97, 0.9, "forced fit", transform=ax.transAxes, horizontalalignment="right", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
            
    # results parameters
        
    ax.text(right*1.01, bottom+(top-bottom)*0.97, results_parameters_3, transform=fig.transFigure, horizontalalignment="left", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
    ax.text(right*1.01, bottom+(top-bottom)*0.97, results_parameters_4, transform=fig.transFigure, horizontalalignment="left", verticalalignment="top", fontsize=sample_plot_details.annotation_fontsize)
          
   
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([np.min(q)/1.1, np.max(q)*1.1])
    ax.set_ylim([np.min(y)/1.1, np.max(y)*1.1])
        
        
    # tick appearances
        
    ax.tick_params(which="major", labelsize=sample_plot_details.tick_fontsize, width=sample_plot_details.tick_width, length=sample_plot_details.tick_length, direction="out")
    ax.tick_params(which="minor", labelsize=sample_plot_details.minor_tick_fontsize, width=sample_plot_details.minor_tick_width, length=sample_plot_details.minor_tick_length, direction="out", labelleft=False)
    ax.xaxis.set_ticks([0.02, 0.04, 0.07, 0.1, 0.2])
    ax.xaxis.set_ticklabels(["0.02", "0.04", "0.07", "0.1", "0.2"])
        
    plt.savefig("Guinier_Porod_fit_NIMROD000%i.png" %(file_start))
        
    plt.close()
        
    return
        


###########     FUNCTION to plot_radius_of_gyration     ##########

def plot_radius_of_gyration(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results):
    
    filename = "%s_%s__Radius_of_Gyration_vs_" %(sample_plot_details.sample, sample_plot_details.action)
    ylabel = "Radius of Gyration ($\AA$)"
    
        
    # color      [color, alpha]
    color = ["r", 0.5]
        
    yscale = [0, 150]   # To be changed by user
        
    Rg1 = sample_guinier_porod_results.Rg
    Rg1_err	= sample_guinier_porod_results.Rg_err
        
    for i in range(sample_guinier_porod_results.number_of_files):
        if Rg1[i]==0:
            Rg1[i] = -100
            Rg1_err[i] = -1
                
    plot = "scatter_errorbar"
    y = [Rg1, Rg1_err]
        
    plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
        
    
    
    return


###########     FUNCTION to plot_s_parameter     ##########

def plot_s_parameter(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results):
    
    filename = "%s_%s__s-parameter_vs_" %(sample_plot_details.sample, sample_plot_details.action)
    ylabel = "s-Parameter"
        
        
    # color      [color, alpha]
    color = [["r", 0.5]]
        
    yscale		= [0,3]   # To be changed by user
        
    s1 = sample_guinier_porod_results.s
    s1_err = sample_guinier_porod_results.s_err
       
        
    for i in range(sample_guinier_porod_results.number_of_files):
        if s1[i]==0:
            s1[i] = -100
            s1_err[i] = -1
                
    plot = "scatter_errorbar"
    y = [s1, s1_err]

    plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
    
    return
    
    


###########     FUNCTION to plot_d_parameter     ##########

def plot_d_parameter(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results):
    
    filename = "%s_%s__d-parameter_vs_" %(sample_plot_details.sample, sample_plot_details.action)
    ylabel = "d-Parameter"
    
        
    # color      [color, alpha]
    color = [["r", 0.5]]
        
    yscale		= [3.4,4.4]   # To be changed by user
        
    d1 = sample_guinier_porod_results.d
    d1_err = sample_guinier_porod_results.d_err
       
        
    for i in range(sample_guinier_porod_results.number_of_files):
        if d1[i]==0:
            d1[i] = -100
            d1_err[i] = -1
                
    plot = "scatter_errorbar"
    y = [d1, d1_err]
        
    plot_python_analysis_results(sample_plot_details, sample_gudrun_results, y, yscale, ylabel, filename, plot, color)
                
    
    return



############################################


'''
##########################################
##                                                                                         ##
##                                  RUNNING FITS                                  ##
##                                                                                         ##
##########################################
'''

	
###########     FUNCTION to check if fit parameters changed     ##########

def number_of_fit_parameters_changed(result_parameters, start_parameters):

	fit_param_changed	= 0
	for p in range(len(result_parameters)):
		if result_parameters[p] != start_parameters[p]:
			fit_param_changed += 1

	print("fit parameters changed: %i" %(fit_param_changed))

	return fit_param_changed


###########     FUNCTION to check if fit uncertainties are real numbers     ##########

def fit_uncertainties_are_real_numbers(result_uncertainties):

	for p in range(len(result_uncertainties)):
		if math.isnan(result_uncertainties[p]) or math.isinf(result_uncertainties[p]):
			print("WARNING\nnan in fit uncertainties")
			return False

	return True


###########     FUNCTION to get fraction of residuals that are outside uncertainties     ##########

def fraction_of_residuals_outside_uncertainties(residuals, uncertainties, factor=1):

	residuals_outside_uncertainties	= 0
	for j in range(len(residuals)):
		if np.abs(residuals[j]) > factor*uncertainties[j]:
			residuals_outside_uncertainties	+= 1

	residuals_outside_uncertainties	= residuals_outside_uncertainties/len(residuals)
	print("residuals > %g x uncertainty: %.0f%%" %(factor, 100. * residuals_outside_uncertainties))

	return residuals_outside_uncertainties


###########     FUNCTION to force_fit     ##########

def find_previous_fit(i, p_0, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files):
    
    s1 = None
    Rg1 = None
    s1_err = None
    Rg1_err = None
    
    # find last scan that was fitted, skipping over exceptions
    
    if reverse_files:
        k = i + 1
        while (k < sample_guinier_porod_results.number_of_files) and (sample_gudrun_results.file_start[k] in sample_plot_details.exceptions):
            k += 1
    else:
        k = i - 1
        while (k >= 0) and (sample_gudrun_results.file_start[k] in sample_plot_details.exceptions):
            k -= 1
            
    if (k < 0) or (k >= sample_guinier_porod_results.number_of_files) or (sample_gudrun_results.file_start[i] in sample_guinier_porod_input.files_manual_start_parameters):
        
        popt = p_0
        perr = [parameter * 0.005 for parameter in popt]
        
    else: 
        
        popt = [sample_guinier_porod_results.G[k], sample_guinier_porod_results.s[k], sample_guinier_porod_results.Rg[k], sample_guinier_porod_results.d[k]]
        perr = [sample_guinier_porod_results.G_err[k], sample_guinier_porod_results.s_err[k], sample_guinier_porod_results.Rg_err[k], sample_guinier_porod_results.d_err[k]]
        
        s1 = sample_guinier_porod_results.s[k]
        Rg1 = sample_guinier_porod_results.Rg[k]
        s1_err = sample_guinier_porod_results.s_err[k]
        Rg1_err = sample_guinier_porod_results.Rg_err[k]

        
    return popt, perr, s1, Rg1, s1_err, Rg1_err


###########     FUNCTION to force_fit     ##########

def force_fit(i, p_0, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files):

	popt, perr, s1, Rg1, s1_err, Rg1_err = find_previous_fit(i, p_0, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files)
	if s1 is not None:
		sample_guinier_porod_results.s[i]		= s1
		sample_guinier_porod_results.Rg[i]		= Rg1
		sample_guinier_porod_results.s_err[i]	= s1_err
		sample_guinier_porod_results.Rg_err[i]	= Rg1_err

	sample_guinier_porod_results.forced_fit[i]	= 1
	print("\nForcing fit\n")

	return popt, perr

	
###########     FUNCTION to run guinier-porod fit     ##########

def run_guinier_porod_fit(i, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results, sample_plot_details, reverse_files=False):

	sample_guinier_porod_results.fit_tag[i]	= 1
    
	print("Guinier-Porod fit")
	print(sample_guinier_porod_temporary_fit_parameters.p0_normal)

	# set fit parameters
	file_start					= sample_gudrun_results.file_start[i]

	x								= sample_guinier_porod_temporary_fit_parameters.q_fit
	y								= np.log(sample_guinier_porod_temporary_fit_parameters.y_fit)
	yerr							= sample_guinier_porod_temporary_fit_parameters.yerr_fit/sample_guinier_porod_temporary_fit_parameters.y_fit
	weight						= yerr/y
	q_min						= sample_guinier_porod_input.q_plot_min
	q_max						= sample_guinier_porod_input.q_plot_max
	
	G, s, Rg, d					= sample_guinier_porod_temporary_fit_parameters.p0_normal
	p_start							= [G, s, Rg, d]

	[G_prev, s_prev, Rg_prev, d_prev], perr_prev, s1_prev, Rg1_prev, s1_err_prev, Rg1_err_prev = find_previous_fit(i, p_start, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files)

	start_residuals	= guinier_porod_function(x, G, s, Rg, d, q_min, q_max) - y
	previous_residuals	= guinier_porod_function(x, G_prev, s_prev, Rg_prev, d_prev, q_min, q_max) - y

	# Check that fit parameters actually change
	count_rerun				= 0
	rerun						= True
	# record all fit results and choose best one, in case reruns need to be attempted


	popt_list				= []
	perr_list				= []
	fit_residuals_list	= []
	while rerun:
	# run fit
		try:
			popt, pcov = curve_fit(lambda x, G, s, Rg, d: guinier_porod_function(x, G, s, Rg, d, q_min, q_max), x, y, p0=[G, s, Rg, d], sigma=weight, maxfev=10000) #absolute_sigma=True
			print(pcov)
			perr 			= np.sqrt(np.diag(pcov))
		except(RuntimeError):
			popt	= [G* 1.0001, s* 1.0001, Rg* 1.0001, d* 1.0001]
			perr	= [G* 0.1, np.max([s* 0.1, 0.3]), Rg* 0.1, d* 0.1]

		# Check fit results
		fit_residuals = np.zeros(len(x)) + 1e10
		print("\n\niteration %i" %(count_rerun + 1))
		if number_of_fit_parameters_changed(popt, [G, s, Rg, d]) > 2:
			if fit_uncertainties_are_real_numbers(perr):
				G_fit, s_fit, Rg_fit, d_fit	= popt
				fit_residuals	= guinier_porod_function(x, G_fit, s_fit, Rg_fit, d_fit, q_min, q_max) - y
				if fraction_of_residuals_outside_uncertainties(fit_residuals, yerr, factor=2.) < 1./3.:
					rerun = False
					if sum_of_squares_larger_for_first_argument(fit_residuals, start_residuals):
						print("use start parameters instead of fit")
						popt, perr = force_fit(i, p_start, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files)
						break

		if rerun:
			popt_list.append(popt)
			perr_list.append(perr)
			fit_residuals_list.append(fit_residuals)

			if count_rerun < 2:
				count_rerun	+= 1
				print("run GP fit again")
				if Rg <= 2./q_max:
					Rg = np.max([1.+2./q_max, 1.+1./q_max*np.sqrt((d-s)*(3-s)/2.)])
					print("new Rg = %.1f" %(Rg))
				else:
					if reverse_files:
						G		= G * 1.05
						s		= s * 1.05	
						Rg	= Rg * 0.95
					else:
						G		= G * 0.95
						s		= s * 0.95	
						Rg	= Rg * 1.05
			else:
				print("stop GP fit")
				# check if any previous fit attempt yielded smaller residuals
				for c in range(count_rerun):
					if sum_of_squares_larger_for_first_argument(fit_residuals, fit_residuals_list[c]):
						popt			 = popt_list[c]
						perr 				= perr_list[c]
						fit_residuals = fit_residuals_list[c]
				# check if fit of preious data set yields smaller residuals for current data set
				if sum_of_squares_larger_for_first_argument(fit_residuals, previous_residuals):
					print("use previous results instead of fit")
					popt, perr = force_fit(i, p_start, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results, reverse_files)
				break

    	
	sample_guinier_porod_results.d[i]						= popt[3]
	sample_guinier_porod_results.G[i]					= popt[0]
	sample_guinier_porod_results.s[i]					= popt[1]
	sample_guinier_porod_results.Rg[i]				= popt[2]
	sample_guinier_porod_results.d_err[i]				= perr[3]
	sample_guinier_porod_results.G_err[i]			= np.min([perr[0], popt[0]*0.2])
	sample_guinier_porod_results.s_err[i]			= np.min([perr[1], 0.3])
	sample_guinier_porod_results.Rg_err[i]			= np.min([perr[2], popt[2]*0.2])
	
    
	sample_guinier_porod_results.porod_range[0][i]	= q_min
	sample_guinier_porod_results.porod_range[1][i]	= q_max

	print("%.2f, %.2f, %.2f, %.2f" %(sample_guinier_porod_results.G[i], sample_guinier_porod_results.s[i], sample_guinier_porod_results.Rg[i], sample_guinier_porod_results.d[i]))


	return popt
	
	


###########     FUNCTION to call_guinier_porod_analysis_loop     ##########

def call_guinier_porod_analysis_loop(i, sample_material_properties, sample_analysis_directories, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_gudrun_results, sample_guinier_porod_results, sample_plot_details, reverse_files=False):
    
    file_start = sample_gudrun_results.file_start[i]
    
    if "sequence" in sample_plot_details.action:

        parameter = sample_gudrun_results.parameter[i]

    else:

        sample = sample_gudrun_results.sample[i]
        
    plot = True
    
    python_fits_read_neutron_data_from_file(i, [sample_guinier_porod_input.q_plot_min, sample_guinier_porod_input.q_plot_max], sample_analysis_directories, sample_guinier_porod_input, sample_gudrun_results, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_results)
    
    if file_start in sample_plot_details.exceptions:
        
        plot = False
        
    elif file_start in sample_guinier_porod_input.files_to_fit_manually:
        
        assign_manual_guinier_porod_fit_parameters(i, sample_gudrun_results, sample_guinier_porod_input, sample_guinier_porod_results)
        
    else:
        
        if file_start in sample_guinier_porod_input.files_manual_start_parameters:
            assign_manual_start_parameters(file_start, sample_guinier_porod_input, sample_guinier_porod_temporary_fit_parameters)
            
        else:
            if reverse_files:
                k = i+1
                while (k < sample_guinier_porod_results.number_of_files) and (sample_gudrun_results.file_start[k] in sample_plot_details.exceptions):
                    k += 1
            else:
                k = i-1
                while (k >= 0) and (sample_gudrun_results.file_start[k] in sample_plot_details.exceptions):
                    k -= 1
            if (k >= 0) and (k < sample_guinier_porod_results.number_of_files):
                # take results from last file's fit, if there is one
                
                sample_guinier_porod_temporary_fit_parameters.p0_porod	= [sample_guinier_porod_results.G[k], sample_guinier_porod_results.d[k]]
                sample_guinier_porod_temporary_fit_parameters.p0_normal= [sample_guinier_porod_results.G[k], sample_guinier_porod_results.s[k], sample_guinier_porod_results.Rg[k], sample_guinier_porod_results.d[k]]
                
                if sample_guinier_porod_temporary_fit_parameters.p0_normal[1] < 0.1:
                    sample_guinier_porod_temporary_fit_parameters.p0_normal[0]	= 1000.
                    sample_guinier_porod_temporary_fit_parameters.p0_normal[1]	= 0.1
   
                    
    # Fit Data
    
    popt  = run_guinier_porod_fit(i, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results, sample_plot_details, reverse_files=reverse_files)
    sample_guinier_porod_temporary_fit_parameters.p0_normal= deepcopy(popt)
    
        
    if sample_guinier_porod_temporary_fit_parameters.p0_normal[1] < 0.1:
        sample_guinier_porod_temporary_fit_parameters.p0_normal[0]	= 1000.
        sample_guinier_porod_temporary_fit_parameters.p0_normal[1]	= 0.1
            
   
    # Plot fit
    if plot:
        plot_individual_guinier_porod_fits(i, sample_plot_details, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
        
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

def save_guinier_porod_results_to_file(input, sample_plot_details, sample_guinier_porod_input, sample_gudrun_results, sample_guinier_porod_results):
    
    if "sequence" in sample_plot_details.action:
        
        np.savetxt("GP_analysis_Porod_range_%s_%s.txt" %(input["sample"], input["action"]), np.transpose([sample_gudrun_results.parameter, sample_guinier_porod_results.porod_range[0], sample_guinier_porod_results.porod_range[1]]), fmt=['%1.4e','%1.4e','%1.4e'], newline=os.linesep, header="Parameter, Porod_range_left, Porod_range_right")
        
        # set up variables to save results to file
        
        Save_Header = ""
        Save_Data = []
        Save_Format = []
        
        # sort results into variables
        
        Save_Header += "parameter, "
        Save_Data.append(sample_gudrun_results.parameter)
        Save_Format.append('%1.4e')
        Save_Header		+= "radius_of_gyration, "
        Save_Data.append(sample_guinier_porod_results.Rg)
        Save_Format.append('%1.4e')
        Save_Header		+= "Rg_err, "
        Save_Data.append(sample_guinier_porod_results.Rg_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "s_parameter, "
        Save_Data.append(sample_guinier_porod_results.s)
        Save_Format.append('%1.4e')
        Save_Header		+= "s_err, "
        Save_Data.append(sample_guinier_porod_results.s_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "d_parameter, "
        Save_Data.append(sample_guinier_porod_results.d)
        Save_Format.append('%1.4e')
        Save_Header		+= "d_err, "
        Save_Data.append(sample_guinier_porod_results.d_err)
        Save_Format.append('%1.4e')
        
        np.savetxt("GP_analysis_%s_%s.txt" %(input["sample"], input["action"]), np.transpose(Save_Data), fmt=Save_Format, header=Save_Header)
        
        
    else:
        
        np.savetxt("GP_analysis_Porod_range_%s_%s.txt" %(input["sample"], input["action"]), np.transpose([sample_gudrun_results.sample, sample_guinier_porod_results.porod_range[0], sample_guinier_porod_results.porod_range[1]]), fmt=['%1.4e','%1.4e','%1.4e'], newline=os.linesep, header="Sample, Porod_range_left, Porod_range_right")
        
        # set up variables to save results to file
        
        Save_Header = ""
        Save_Data = []
        Save_Format = []
        
        # sort results into variables
        
        Save_Header += "sample, "
        Save_Data.append(sample_gudrun_results.sample)
        Save_Format.append('%1.4e')
        Save_Header		+= "radius_of_gyration, "
        Save_Data.append(sample_guinier_porod_results.Rg)
        Save_Format.append('%1.4e')
        Save_Header		+= "Rg_err, "
        Save_Data.append(sample_guinier_porod_results.Rg_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "s_parameter, "
        Save_Data.append(sample_guinier_porod_results.s)
        Save_Format.append('%1.4e')
        Save_Header		+= "s_err, "
        Save_Data.append(sample_guinier_porod_results.s_err)
        Save_Format.append('%1.4e')
        Save_Header		+= "d_parameter, "
        Save_Data.append(sample_guinier_porod_results.d)
        Save_Format.append('%1.4e')
        Save_Header		+= "d_err, "
        Save_Data.append(sample_guinier_porod_results.d_err)
        Save_Format.append('%1.4e')
        
        np.savetxt("GP_analysis_%s_%s.txt" %(input["sample"], input["action"]), np.transpose(Save_Data), fmt=Save_Format, header=Save_Header)
        
    return


############################################

'''
##########################################
##                                                                                         ##
##                                 MAIN FUNCTION                                ##
##                                                                                         ##
##########################################
'''


###########     FUNCTION to run_guinier_porod_analysis     ##########

def run_guinier_porod_analysis(input, sample_material_properties, sample_analysis_directories, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input):
	print("\n\n-------------------\n\nRunning Guinier-Porod analysis.\n\n")

	# Initialise arrays for fit results
	try:
		n	= len(sample_gudrun_results.file_start)
	except(TypeError):
		n	= 1
	sample_guinier_porod_results = guinier_porod_results(n)

	
	loop_range	= range(sample_guinier_porod_results.number_of_files-1, -1, -1)
	reverse_files=True
	
	for i in loop_range:
		# Initialise temporary storage for data and start parameters
		sample_guinier_porod_temporary_fit_parameters = guinier_porod_temporary_fit_parameters()
		# Call fit	
		call_guinier_porod_analysis_loop(i, sample_material_properties, sample_analysis_directories, sample_guinier_porod_temporary_fit_parameters, sample_guinier_porod_input, sample_gudrun_results, sample_guinier_porod_results, sample_plot_details, reverse_files=reverse_files)
	
	print("\n\n------------------------------------\n\nGuinier Porod Summary\n")
	print("%i Python fits\n" %(np.count_nonzero(sample_guinier_porod_results.fit_tag==1)))
	for i in range(sample_guinier_porod_results.number_of_files):
		if sample_guinier_porod_results.forced_fit[i] == 1:
			print(sample_gudrun_results.file_start[i], "\tforced fit")
		if sample_guinier_porod_results.fit_tag[i] == 0:
			print(sample_gudrun_results.file_start[i], "\tmanual fit")
		if sample_guinier_porod_results.fit_tag[i] == 2:
			print(sample_gudrun_results.file_start[i], "\tnot analysed")
	print("\n\n")
	
	# Save results to file
	save_guinier_porod_results_to_file(input, sample_plot_details, sample_guinier_porod_input, sample_gudrun_results, sample_guinier_porod_results)
	
	# Plot results
	plot_radius_of_gyration(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
	plot_s_parameter(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
	plot_d_parameter(sample_plot_details, sample_guinier_porod_input, sample_guinier_porod_results, sample_gudrun_results)
