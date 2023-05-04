# Model variables

vasculature_time        = 1
th                      = 1e-3
tha                     = th
xh                      = 5e-3
xha                     = xh
yh                      = 1e-4
dM                      = 1e-4
dE                      = 5e-5
phiM                    = 5e-4
phiE                    = 5e-4
dmmp                    = 1e-3
theta                   = 0.195
Lambda                  = 0.1
gamma1                  = 1
gamma2                  = 1
doublingTimeE           = 2000 # PDF: 2000
doublingTimeM           = 3000 # PDF: 3000
gridsize                = 10
totalTime               = 24002 # must contain step 0 and final step of batch. (I think this is never being used)
patchsize               = 3
middlePoint             = round(gridsize/2)
E1                      = 0.5461
E2                      = 0.2553
E3                      = 0.1986
single_cell_survival    = 1         # 5e-4 #probability of survival for single cells in the vasculature
cluster_survival        = 1 #2.5e-2    # probability of survival for clusters in the vasculature
dissagreggation_prob    = 0.5       # probability of dissagreggation in vasculature
carrying_capacity       = 4

# These values bellow change according to the simulation we want to recreate
mesenchymal_proportion = 0.6 # PDF: 0.6
epithelial_proportion = 0.4 # PDF: 0.4
n_center_points_for_tumor = 97 #PDF: 97
n_center_points_for_Vessels = 200 #PDF 200 # vessels will not be in this center points

# Variable to stop the saving of ecm and mmp2 if it is NOT a Batch run analysis
isBatchRun = True

# Name of the directorie to save the data
parent_dir = "Simulations"
imagesFolder_utils = "Visual analysis"
gridsize_utils = 201

