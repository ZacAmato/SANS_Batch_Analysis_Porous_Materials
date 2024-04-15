import os
import numpy as np
from CLASSES import *
from plot_settings import *


###########     FUNCTION to read data from file     ##########

def python_fits_read_neutron_data_from_file(i, fit_range, sample_analysis_directories, sample_input, sample_gudrun_results, sample_plot_details, sample_temporary_fit_parameters, sample_results):
	# Exception handling required to catch cases with only one data point, where indexing won't work.
	file_start	= sample_gudrun_results.file_start[i]
	
		
	if file_start in sample_plot_details.exceptions:
		print("---------------------------------------\nFile:%i\t%i of %i\not plotted" %(file_start,(i+1), sample_results.number_of_files))
		try:
			sample_results.fit_tag[i]		= 2
		except(AttributeError):
			pass
	else:
		print("---------------------------------------\nFile:%i\t%i of %i" %(file_start, (i+1), sample_results.number_of_files))
		
		#Get values
		os.chdir(sample_analysis_directories.Gudrun_directory)
		q, y, yerr	= np.loadtxt("NIMROD000%i.mint01" %(file_start), unpack=True)
		os.chdir(sample_analysis_directories.Python_directory)

		#Slice arrays.
		r		= np.where((q >= sample_input.q_plot_min) & (q <= sample_input.q_plot_max) & (y > 0))
		r_fit	= np.where((q >= fit_range[0]) & (q <= fit_range[1]) & (y > 0))

		sample_temporary_fit_parameters.q	 		= q[r]
		sample_temporary_fit_parameters.y 			= y[r]
		sample_temporary_fit_parameters.yerr		= yerr[r]
		sample_temporary_fit_parameters.q_fit		= q[r_fit]
		sample_temporary_fit_parameters.y_fit 		= y[r_fit]
		sample_temporary_fit_parameters.yerr_fit	= yerr[r_fit]
		
	return
	


###########     FUNCTION to smooth_data_over_window     ##########

def smooth_data_over_window(data, window_width):
	
	# Note: window_width must be odd integer
	
	#~ print("\nSmoothing data")
	
	data_smoothed	= np.zeros(np.size(data))
	window_index	= int((window_width-1)/2)
	
	for k in range(np.size(data)):
		start = k - window_index
		stop = k + window_index + 1

		if(start < 0):
			stop += start
			start = 0

		if(stop > np.size(data)):
			start += stop - np.size(data)
			stop = np.size(data)

		data_smoothed[k]	= np.average(data[start:stop])

	return data_smoothed

'''
			# Compute rate of change.

'''
###########     FUNCTION to get_rate_of_change     ##########

def get_rate_of_change(data):
	
	rate_of_change 		= np.zeros(np.size(data))
	for k in range(np.size(data)):
		rate_of_change[k] = data[k] - data[k-1]

	return rate_of_change


###########     FUNCTION to check if sum_of_squares_larger_for_first_argument     ##########

def sum_of_squares_larger_for_first_argument(array_1, array_2):

	if np.sum(np.square(array_1)) > np.sum(np.square(array_2)):
		return True
	else:
		return False


###########     FUNCTION to read_gudrun_analysis_input_file     ##########

def read_gudrun_analysis_input_file(input_file="Input_Parameters_Python_Analysis.csv"):
	
	# Read input parameters for analysis from file
	sample_input_parameters	= np.loadtxt(input_file, dtype="str", delimiter=",")
	# Sort input parameters into dictionary
	input			= {}
	for line in sample_input_parameters:
		input[line[0]]	= line[1]

	# Activate full print-outs?		(True or False)
	verbose = False

	# Get directories
	if input["work_directory"] == ".":
		input["work_directory"]		= os.getcwd()
	if "." in input["results_directory"]:
		input["results_directory"]	= input["results_directory"].replace(".", os.getcwd())
	
	return input, verbose
	
	
###########     FUNCTION to update input dictionary    ##########

def update_input_dictionary(parameter, keys, dictionary, data_source, data_type=None):
	if data_source[2] == parameter:
		for i in range(len(keys)):
			if data_source[3+i] != "":
				key	= keys[i]
				data	= data_source[3+i]
				if data_type == "int":
					data = int(data)
				elif data_type == "float":
					data = float(data)
				else:
					print("\nunknown data type specified\n")
				
				dictionary[key]	= data
		
	return dictionary	
	
	
###########     FUNCTIONs to set python analysis input    ##########

def set_guinier_porod_input(sample_guinier_porod_input, line):
    
    if line[2] == "q plot min" and line[3] != "":
        sample_guinier_porod_input.q_plot_min = float(line[3])
    elif line[2] == "q plot max" and line[3] != "":
        sample_guinier_porod_input.q_plot_max = float(line[3])
        
        
    # Manual fits    
    elif line[2] == "files to fit manually" and line[3] != "":
        for file in line[3:]:
            try:
                sample_guinier_porod_input.files_to_fit_manually.append(int(file))
            except(ValueError):
                pass
    sample_guinier_porod_input.manual_fit_d = update_input_dictionary(parameter="manual fit d", keys=sample_guinier_porod_input.files_to_fit_manually, dictionary=sample_guinier_porod_input.manual_fit_d, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_fit_G_low_Q = update_input_dictionary(parameter="manual fit G low Q", keys=sample_guinier_porod_input.files_to_fit_manually, dictionary=sample_guinier_porod_input.manual_fit_G_low_Q, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_fit_s_low_Q = update_input_dictionary(parameter="manual fit s low Q", keys=sample_guinier_porod_input.files_to_fit_manually, dictionary=sample_guinier_porod_input.manual_fit_s_low_Q, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_fit_Rg_low_Q = update_input_dictionary(parameter="manual fit Rg low Q", keys=sample_guinier_porod_input.files_to_fit_manually, dictionary=sample_guinier_porod_input.manual_fit_Rg_low_Q, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_fit_uncertainty = update_input_dictionary(parameter="manual fit uncertainty in %", keys=sample_guinier_porod_input.files_to_fit_manually, dictionary=sample_guinier_porod_input.manual_fit_uncertainty, data_source=line, data_type="float")
    

    # Manual start parameters
    if line[2] == "files to set fit start parameters manually" and line[3] != "":
        for file in line[3:]:
            try:
                sample_guinier_porod_input.files_manual_start_parameters.append(int(file))
            except(ValueError):
                pass
    sample_guinier_porod_input.manual_start_d = update_input_dictionary(parameter="manual start d", keys=sample_guinier_porod_input.files_manual_start_parameters, dictionary=sample_guinier_porod_input.manual_start_d, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_start_G_low_Q = update_input_dictionary(parameter="manual start G low Q", keys=sample_guinier_porod_input.files_manual_start_parameters, dictionary=sample_guinier_porod_input.manual_start_G_low_Q, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_start_s_low_Q = update_input_dictionary(parameter="manual start s low Q", keys=sample_guinier_porod_input.files_manual_start_parameters, dictionary=sample_guinier_porod_input.manual_start_s_low_Q, data_source=line, data_type="float")
    sample_guinier_porod_input.manual_start_Rg_low_Q = update_input_dictionary(parameter="manual start Rg low Q", keys=sample_guinier_porod_input.files_manual_start_parameters, dictionary=sample_guinier_porod_input.manual_start_Rg_low_Q, data_source=line, data_type="float")
    
    return sample_guinier_porod_input
        
    

def set_porod_constant_input(sample_porod_constant_input, line):
	if line[2] == "files to set q fit range manually" and line[3] != "":
		for file in line[3:]:
			try:
				sample_porod_constant_input.files_manual_fit_range.append(int(file))
			except(ValueError):
				pass
	sample_porod_constant_input.q_fit_range_min = update_input_dictionary(parameter="q fit range min", keys=sample_porod_constant_input.files_manual_fit_range, dictionary=sample_porod_constant_input.q_fit_range_min, data_source=line, data_type="float")
	sample_porod_constant_input.q_fit_range_max = update_input_dictionary(parameter="q fit range max", keys=sample_porod_constant_input.files_manual_fit_range, dictionary=sample_porod_constant_input.q_fit_range_max, data_source=line, data_type="float")
	if line[2] == "number of points to fit over" and line[3] != "":
		sample_porod_constant_input.number_of_points_to_fit_over = int(line[3])
	elif line[2] == "q plot min" and line[3] != "":
		sample_porod_constant_input.q_plot_min = float(line[3])
	elif line[2] == "q plot max" and line[3] != "":
		sample_porod_constant_input.q_plot_max = float(line[3])
	return sample_porod_constant_input
	
	
	
###########     FUNCTION to write input parameters to file for checking     ##########

def write_input_checking_file(sample_porod_constant_input, sample_guinier_porod_input):
	output = []
	
	output.append("\nGuinier-Porod\n\n")
	output.append("q plot min\t%s\n" %(	sample_guinier_porod_input.q_plot_min))
	output.append("q plot max\t%s\n" %(sample_guinier_porod_input.q_plot_max))
	output.append("files to set fit start parameters manually\t%s\n" %(sample_guinier_porod_input.files_manual_start_parameters))
	output.append("manual start G low Q\t%s\n" %(sample_guinier_porod_input.manual_start_G_low_Q))
	output.append("manual start s low Q\t%s\n" %(sample_guinier_porod_input.manual_start_s_low_Q))
	output.append("manual start Rg low Q\t%s\n" %(sample_guinier_porod_input.manual_start_Rg_low_Q))	

	output.append("\nPorod Constant\n\n")
	output.append("files to set q fit range manually\t%s\n" %(sample_porod_constant_input.files_manual_fit_range))
	output.append("q fit range min\t%s\n" %(sample_porod_constant_input.q_fit_range_min))
	output.append("q fit range max\t%s\n" %(sample_porod_constant_input.q_fit_range_max))
	output.append("number of points to fit over\t%s\n" %(sample_porod_constant_input.number_of_points_to_fit_over))
	output.append("q plot min\t%s\n" %(sample_porod_constant_input.q_plot_min))
	output.append("q plot max\t%s\n" %(sample_porod_constant_input.q_plot_max))
	
	
	f	= open("check_input_file.txt", 'w')
	f.writelines(output)
	f.close()

	return
		
		
###########     FUNCTION to read_python_analysis_input_file     ##########

def read_python_analysis_input_file(input_file="Input_Parameters_Python_Analysis.csv"):
	
	# Read input parameters for analysis from file
	python_input_parameters			= np.loadtxt(input_file, dtype="str", delimiter=",")

	# Initialise input parameter classes
	sample_porod_constant_input = porod_constant_input()
	sample_guinier_porod_input		= guinier_porod_input()

	
	
	for line in python_input_parameters:
		# Exception handling to udpdate only those input parameters that are specified in input file.
		# All others are left at default values.
		try:
			# Guinier-Porod input
			if line[1] == "Guinier-Porod":
				sample_guinier_porod_input		= set_guinier_porod_input(sample_guinier_porod_input, line)

			elif line[1] == "Porod Constant":
				sample_porod_constant_input = set_porod_constant_input(sample_porod_constant_input, line)
			
		
		except (IndexError):
			pass
			
	# Write all input parameters to output file for checking.
	write_input_checking_file(sample_porod_constant_input, sample_guinier_porod_input)

	return 	sample_porod_constant_input, sample_guinier_porod_input


###########     FUNCTION to set default python analysis parameters and info    ##########

def set_python_analysis_parameter_and_info(input):
	print("\n" + input["sample"] + "\n\t" + input["action"] + "\n")

	# Set directories & material properties
	sample_analysis_directories							= analysis_directories()
	sample_material_properties							= material_properties()
	
	# Get plot details and gudrun results

	sample_plot_details, sample_gudrun_results	= get_plot_info(input)
	sample_plot_details										= define_color_scheme(sample_plot_details, sample_gudrun_results)
	
	return sample_analysis_directories, sample_material_properties, sample_plot_details, sample_gudrun_results 
