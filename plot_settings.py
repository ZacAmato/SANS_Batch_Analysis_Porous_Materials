import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib as mpl
from matplotlib.colors import *
from CLASSES import *


###########     FUNCTION to get input info and tags for plotting     ##########

def get_plot_info(input):
    
	
    sample_plot_details	= plot_details(
								sample					= input["sample"]					,
								action					= input["action"]					,
								work_directory		= input["work_directory"]		,
								results_directory	= input["results_directory"]	,
								)

	# get exceptions
    try:
    	for file in input["exceptions"].split(","):
    		sample_plot_details.exceptions.append(int(file))
    except (KeyError):
    	pass

	
	# Get files and labels:
        
    if "sequence" in sample_plot_details.action:
            
            file_start, parameter = np.loadtxt("results_%s_%s.txt" %(sample_plot_details.sample, sample_plot_details.action),  unpack=1)
            
            try:
                if len(file_start) > 1:
                    pass
                sample_gudrun_results	= gudrun_results(
                                   file_start			= file_start		,
											parameter				= parameter				,
											)
            except(TypeError):
                sample_gudrun_results	= gudrun_results(
										file_start			= [file_start]	,
										parameter				= [parameter]				,
										)
                
            return sample_plot_details, sample_gudrun_results
            
    else:
        
        file_start, sample = np.loadtxt("results_%s_%s.txt" %(sample_plot_details.sample, sample_plot_details.action),  unpack=1)
        
        try:
                if len(file_start) > 1:
                    pass
                sample_gudrun_results	= gudrun_results2(
                                   file_start			= file_start		,
											sample				= sample				,
											)
        except(TypeError):
            sample_gudrun_results	= gudrun_results2(
									file_start			= [file_start]	,
									sample				= [sample]				,
									)
                
        return sample_plot_details, sample_gudrun_results
            


###########     FUNCTION to set the colour scheme     ##########

def define_color_scheme(sample_plot_details, sample_gudrun_results):
	
    # Set colors, labels & ranges, depending on action to be plotted
    if "sequence" in sample_plot_details.action:
        
            sample_plot_details.color_bar_label = "Parameter"                 
            sample_plot_details.color_range			= [80., 180.]  # To be changed by user
            tick_number										= 11   # To be changed by user
            sample_plot_details.color_values	= sample_gudrun_results.parameter
            

            sample_plot_details.color_norm 			= mpl.colors.Normalize(
																		vmin=sample_plot_details.color_range[0]	, 
																		vmax=sample_plot_details.color_range[-1]
																		)
            sample_plot_details.color_bounds 		= np.linspace(
																sample_plot_details.color_range[0]				, 
																sample_plot_details.color_range[-1]			, 
																500
																)
            sample_plot_details.color_plot_ticks	= np.linspace(
																sample_plot_details.color_range[0]				, 
																sample_plot_details.color_range[-1]			, 
																tick_number
																)

            return sample_plot_details

    else:

    	return sample_plot_details
    	  
