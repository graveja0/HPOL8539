* ABOUT THIS CODE COLLECTION
This collection of code replicates the simulations, figures, and examples (except for the real-data example) in the paper.  If any aspect of the paper results (other than the real data application) is missing, please contact the authors so this can be rectified.

Bounding, An Accessible Method for Estimating Principal Causal Effects, Examined and Explained
by Luke Miratrix, Jane Furey, Avi Feller, Todd Grindal & Lindsay C. Page

There is also code demonstrating how to use the functions on other data.

For questions or comments (or found errors!) please contact
Luke Miratrix at lmiratrix <at> gse.harvard.edu


LIST OF FILES
ECHS Project.Rproj
    RStudio project file
    
README.txt
    This document
    
bound_function_library.R
    The functions to calculate bounds, both stratified and not.
    
data_generators.R
    Code for the data generation for the simulation studies
    
demo_bound_package.R
    Illustrating script of the code
    
generate_demo_of_strata_utility_plot.R
    Look at impact of stratification--single example.
    
results
    Scripts will store results here.
    
set_paper_parameters.R
    This script sets up a DGP with specific parameters calibrated to specific scenarios.  Change the flags in this to change the automatically shared contexts of the simulations.
    
study_impact_of_covariate_predict_plot.R
    Main simulation in the paper

study_uncertainty_vs_number_strata.R
    Simulation to examine how uncertainty depends on number of strata.
        
uncertainty_by_strata_plot.R
    Illustrative plot to understand how bounding works.

