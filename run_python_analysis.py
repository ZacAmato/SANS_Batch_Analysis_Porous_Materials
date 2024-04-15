# import my own functions
from analysis_input_file_handling import *
from Guinier_Porod import *
from Porod_Constant import *



# Read input parameters for analysis from file
input, verbose = read_gudrun_analysis_input_file()
sample_porod_constant_input, sample_guinier_porod_input= read_python_analysis_input_file()

# Read and set general sample info
sample_analysis_directories, sample_material_properties, sample_plot_details, sample_gudrun_results = set_python_analysis_parameter_and_info(input)

print(sample_plot_details)

# GUINIER-POROD analysis
run_guinier_porod_analysis(input, sample_material_properties, sample_analysis_directories, sample_plot_details, sample_gudrun_results, sample_guinier_porod_input)

# POROD CONSTANT analysis
run_porod_constant_analysis(input, sample_material_properties, sample_analysis_directories, sample_plot_details, sample_gudrun_results, sample_porod_constant_input)

# END
print("\n\n****************\n\nAnalysis Completed\n\n****************\n\n")
