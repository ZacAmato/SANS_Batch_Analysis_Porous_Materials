'''
ASW data analysis

CLASSES for passing variables between functions

'''

import numpy as np
import os
from datetime import datetime
from analysis_input_file_handling import *
from plot_python_analysis_results import *

###########     Gudrun     ##########


class analysis_parameters(object):
	def __init__(self, work_directory, file_start, sample_counter=0):
		# work directory for python scripts and outputs
		self.work_directory				= work_directory
		# lists of files for Gudrun to analyse (assuming consecutive files belong to one sample)
		self.file_start					= file_start		# first scans on each sample				
		# counters for total number of samples to be analysed and current sample
		self.number_of_samples	= len(file_start)
		self.sample_counter		= sample_counter
		
		

class analysis_results(object):
	def __init__(self, analysed_files=[], analysis_issues=[], analysis_counter=0, issue_counter=0):
		# information to be saved to files to log progress and results
		self.analysed_files			= analysed_files
		self.analysis_issues		= analysis_issues
		# counters
		self.analysis_counter		= analysis_counter
		self.issue_counter			= issue_counter


###########     Gudrun results and plotting     ##########


            
class gudrun_results(object):
                
    def __init__(self,  file_start, parameter):
        self.file_start				= file_start
        self.parameter					= parameter
                             

class gudrun_results2(object):
                
    def __init__(self,  file_start, sample):
        self.file_start				= file_start
        self.sample					= sample
                    
        


class plot_details(object):
	def __init__(self, sample, action, work_directory, results_directory, exceptions=[], plot_title="", color_values=None, color_map=None, color_scale=None, color_range=None, color_norm=None, color_bounds=None, color_plot_ticks=None, color_bar_label=None):
		self.sample				= sample
		self.action					= action
		self.work_directory		= work_directory
		self.results_directory	= results_directory
		self.exceptions			= exceptions
		self.plot_title				= plot_title

		self.color_values		= color_values
		self.color_map			= color_map
		self.color_scale			= color_scale
		self.color_range			= color_range
		self.color_norm			= color_norm
		self.color_bounds		= color_bounds
		self.color_plot_ticks	= color_plot_ticks
		self.color_bar_label	= color_bar_label

		self.label_fontsize			= 20
		self.annotation_fontsize	= 14
		self.tick_fontsize			= 18
		self.tick_width				= 2
		self.tick_length				= 6
		self.minor_tick_fontsize	= 16
		self.minor_tick_width		= 2
		self.minor_tick_length	= 3


###########     Python analysis input file handling     ##########

class analysis_directories(object):
	def __init__(self, Gudrun_directory="Gudrun_Results", Python_directory="python_analysis"):
		self.Gudrun_directory	= os.getcwd() + "/" + Gudrun_directory
		self.Python_directory	= os.getcwd() + "/" + Python_directory

class material_properties(object):
	def __init__(self, atomic_number_density=0.094, scattering_length_density_difference=6.38e-6):              # To be changed by user
		self.atomic_number_density					= atomic_number_density						# (atoms/angstrom^3)
		self.scattering_length_density_difference	= scattering_length_density_difference	# (angstrom^-2)

	
###########     Python analysis parameter and results handling     ##########

# GUINIER-POROD

class guinier_porod_input(object):  
    def __init__(self, q_plot_min=0.012, q_plot_max=0.3, files_to_fit_manually=[], manual_fit_d={}, manual_fit_G_low_Q={}, manual_fit_s_low_Q={}, manual_fit_Rg_low_Q={}, manual_fit_uncertainty={}, files_manual_start_parameters=[], manual_start_d={}, manual_start_G_low_Q={}, manual_start_s_low_Q={}, manual_start_Rg_low_Q={}):
        self.q_plot_min     = q_plot_min
        self.q_plot_max     = q_plot_max
        self.files_to_fit_manually   = files_to_fit_manually
        self.manual_fit_d   = manual_fit_d
        self.manual_fit_G_low_Q     = manual_fit_G_low_Q
        self.manual_fit_s_low_Q   = manual_fit_s_low_Q
        self.manual_fit_Rg_low_Q    = manual_fit_Rg_low_Q
        self.manual_fit_uncertainty     = manual_fit_uncertainty
        self.files_manual_start_parameters  = files_manual_start_parameters
        self.manual_start_d     = manual_start_d
        self.manual_start_G_low_Q   = manual_start_G_low_Q
        self. manual_start_s_low_Q  = manual_start_s_low_Q
        self.manual_start_Rg_low_Q  = manual_start_Rg_low_Q
        
  

class guinier_porod_results(object):
	def __init__(self, number_of_files):
		self.number_of_files	= number_of_files
		self.d					= np.array(np.zeros(number_of_files))
		self.d_err				= np.array(np.zeros(number_of_files))
		self.G					= np.array(np.zeros(number_of_files))
		self.G_err				= np.array(np.zeros(number_of_files))
		self.Rg					= np.array(np.zeros(number_of_files))
		self.Rg_err			= np.array(np.zeros(number_of_files))
		self.s						= np.array(np.zeros(number_of_files))
		self.s_err				= np.array(np.zeros(number_of_files))
		self.fit_tag				= np.array(np.zeros(number_of_files))
		self.forced_fit			= np.array(np.zeros(number_of_files))
		self.porod_range	= [np.zeros(number_of_files), np.zeros(number_of_files)]

class guinier_porod_temporary_fit_parameters(object):
	def __init__(self, q=[], y=[], yerr=[], q_fit=[], y_fit=[], yerr_fit=[], p0_normal=[2000., 1.4, 75., 3.8], p0_porod=[0.0002, 3.77]):
		self.q	 			= q
		self.y 				= y
		self.yerr 			= yerr
		self.q_fit 			= q_fit
		self.y_fit 			= y_fit
		self.yerr_fit		= yerr_fit
		self.p0_normal	= p0_normal			# start values for fit: G, s, Rg, d
		self.p0_porod		= p0_porod			# start values for fit: A, d
		self.fit_attempts	= 0
		self.d					= []
		self.G					= []
		self.Rg					= []
		self.s						= []


# POROD CONSTANT

class porod_constant_input(object):
	def __init__(self, files_manual_fit_range=[], q_fit_range_min={}, q_fit_range_max={}, number_of_points_to_fit_over=15, expand_fit_range=0.2, q_plot_min=0.012, q_plot_max=0.4):
		self.files_manual_fit_range				= files_manual_fit_range
		self.q_fit_range_min						= q_fit_range_min
		self.q_fit_range_max						= q_fit_range_max
		self.number_of_points_to_fit_over	= number_of_points_to_fit_over
		self.expand_fit_range						= expand_fit_range
		self.q_plot_min								= q_plot_min
		self.q_plot_max								= q_plot_max

class porod_constant_results(object):
	def __init__(self, number_of_files):
		self.number_of_files	= number_of_files
		self.SSA					= np.array(np.zeros(number_of_files))
		self.SSA_err				= np.array(np.zeros(number_of_files))
		self.K						= np.array(np.zeros(number_of_files))
		self.K_err					= np.array(np.zeros(number_of_files))+1
		self.q_fit_start			= np.array(np.zeros(number_of_files))
		self.q_fit_end				= np.array(np.zeros(number_of_files))
		self.q_data				= []
		self.y_data				= []
		self.y_data_err			= []
		self.K_line					= []
		self.fit_deviations		= []

class porod_constant_temporary_fit_parameters(object):
	def __init__(self, p0=[], q=[], y=[], yerr=[], q_fit=[], y_fit=[], yerr_fit=[]):
		# self.p0	 			= p0
		#~ self.p1	 			= p1
		self.q	 					= q
		self.y 						= y
		self.yerr		 			= yerr
		self.q_fit 					= q_fit
		self.y_fit 					= y_fit
		self.yerr_fit				= yerr_fit
		self.fit_attempts		= 0
		self.q_fit_ranges		= []
		self.K						= []
		self.K_err					= []
		self.y_err					= []









