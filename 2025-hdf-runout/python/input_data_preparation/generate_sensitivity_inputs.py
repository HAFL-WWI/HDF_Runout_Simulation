######################################################################
# Copyright (C) 2025 BFH
#
# Script for preparing the parameter combinations used in the 
# simulations for the sensitivity analysis .
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

import numpy as np
import pandas as pd
import os
import sys

from SALib.sample.sobol import sample

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_4/code/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

#
# Sample parameter combinations for simulations with forest with variable friction
#

output_base_path = "E:/GIS_Projekte/Paper_4/data/sensitivity_forest"
utilities.ensure_dir(output_base_path)

#Create a problem dictionary. Here we supply the number of variables, the names of each variable, and the bounds of the variables.
problem = {
    "num_vars": 6,
    "names": ["thickness", "density", "height", "type", "mu", "deposit"],
    "bounds": [[0.2, 2], # slide thickness
               [1000, 2000], # moving mass density
               [0.5, 2], #maximum flow height
               [0,3], # slide type
               [0,1], # forest mu
               [0,1] # forest deposit
               ]
}
# Generate parmeter values using the sobol.sample function
param_values = sample(problem, 1024)

param_values[:,3] = np.round(param_values[:,3],0)
param_values[:,4] = np.round(param_values[:,4],0)
param_values[:,5] = np.round(param_values[:,5],0)

param_out_path = os.path.join(output_base_path,"sensitivity_param_values_variable.txt")
np.savetxt(param_out_path, param_values)


#
# Sample parameter combinations for simulations on constant gradient slope with variable friction
#

output_base_path = "E:/GIS_Projekte/Paper_4/data/sensitivity_constant_slope"
utilities.ensure_dir(output_base_path)

problem = {
    "num_vars": 10,
    "names": ["thickness", "density", "height", "type", "friction", "deposit"
             ,"size","slope","dbh","stemdensity"],
    "bounds": [[0.2, 2], # slide thickness
               [1000, 2000], # moving mass density
               [0.5, 2], #maximum flow height
               [0,3], # slide type
               [0,1], # forest friction
               [0,1], # forest deposit
               [0,3], # start size class
               [10,60], # slope gradient 
               [20,50], # DBH
               [200,600] # stem density per ha
               ]
}
# Generate parmeter values using the sobol.sample function
param_values = sample(problem, 1024)

param_values[:,3] = np.round(param_values[:,3],0) # slide type 0/1/2/3
param_values[:,4] = np.round(param_values[:,4],0) # friction 0/1
param_values[:,5] = np.round(param_values[:,5],0) # deposit 0/1
param_values[:,6] = np.round(param_values[:,6],0) # start size class 0:100/1:300/2:600/3:900

param_out_path = os.path.join(output_base_path,"sensitivity_param_values_constant_slope.txt")
np.savetxt(param_out_path, param_values)
