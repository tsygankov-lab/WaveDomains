## WaveDomains
#2D model and analysis of wave domains is cell cortex

run_RD_system.m - main script for running reaction-diffusion (RD) model. Requires specification of parameters, input folder, folder with supplemental functions and output folder

HE_system_2D_noise.m - function for running the 2D version RD model with noise

f.h, h.m - supplemental functions for computing reaction terms in RD model

laplacian.m - supplemental functions for computing laplacian with no-flux boundary conditions for arbitrary shaped cell mask

extract_oocyte_cell_mask.m - script for extracting cell mask from experimental data

compute_textural_features.m - function for computing textural features

wave_direction.m - function for computing wave direction

coherent_distance_map.m - function for computing coherence distance map

WS_segmentation.m - function for watershed segmentation

merge_domains.m - function for domains merging after segmentation

run_2D_square_FEM.m - script for running RD model with Finite Elements Method
