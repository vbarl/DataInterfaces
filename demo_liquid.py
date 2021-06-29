#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------------------------------
Created on Mon May 17 09:06:23 2021
-----------------------------------------------------------------
@purpose : Script for loading different habits from the ARTS DDB
@author  : Vasileios Barlakas
@email   : vasileios.barlakas@gmail.com
           vasileios.barlakas@chalmers.se
-----------------------------------------------------------------
"""
#################################################################
#################################################################
#
# Load libraries
#
import utils,assp
import numpy   as np
#################################################################
#################################################################
#
# Definitions
#
#################################################################
#################################################################
#
# Give database path + initialisation
SSDBpath = '/home/barlakas/WORKAREA/DataBase/SSD'
utils.ssdb_init( SSDBpath )

# Choose orientation
orient = 'totally_random' 


# Choose habit
habID = 25 # liquid spheres 

# Choose sizes
minD = 0. 	    # minimum size to extract
maxD = float("inf") # maximum size to extract

# Choose dmax/dveq selection
sizeparam = 'dmax' # size parameter to apply extraction limits on.
                   # available (param [unit]): 'dmax' [m], 'dveq' [m], 'mass' [kg]

# Choose frequency range    
fmin = 0.           # minimum frequency to extract [Hz]
fmax = float("inf") # maximum frequency to extract
    
# Choose temperature range
tmin = 269. # minimum temperature to extract [K]
tmax = 271. # maximum temperature to extract
    
Si,Mi = assp.assp_import_ssdb( habID, orient, allow_nodata=False,
                         size_range=[minD, maxD], size_type=sizeparam,
                         freq_range=[fmin, fmax],
                         temp_range=[tmin, tmax] )

# The following is what you tried out. For me, it worked just fine. You could add the additional 
# input, i.e. size_range, freq_range, temp_range, exactly as assp.assp_import_ssdb
#
# data = utils.ssdb_import_habit(habID, orient)
