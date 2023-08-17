# BEM-thin-liquid-film-code-data
Accompanying Matlab and c code for "A positivity-preserving numerical method for a thin liquid film on a vertical cylindrical fiber".

B. Kim, H. Ji, A. L. Bertozzi, A. Sadeghpour, and Y. S. Ju. *A positivity-
preserving numerical method for a thin liquid film on a vertical cylindrical fiber*
Under Review in Journal of Computational Physics (2023)

## Figure 2 data and generating code
Data and code for generating Figure 2 is stored in "Figure2". Running "plot_GM_BEM_Fig2_dimensionless.m" will use "Figure2(a)_GM_coarse_data.dat", "Figure2(b)_BEM_coarse_data.dat", and "parameters_used_for_Fig2.dat" to generate Figure 2(a) and 2(b).

## Figure 3 data and generating code
Data and code for generating Figure 3 is stored in "Figure3". Running "plot_GM_BEM_Fig3_dimensionless.m" will use "BEM_coarse_for_Fig3.dat", "GM_coarse_for_Fig3.dat", "GM_fine_for_Fig3.dat", "params_for_coarse_Fig3.dat", and "params_for_fine_Fig3.dat" to generate Figure 3(a) and 3(b).

## Figure 4 data and generating code
Data and code for generating Figure 4 is stored in "Figure4". Running "create_segmented_experimental_image_RP.m" will preprocess use "Experimental_data_RP.png" to produce segmented image. "GM_coarse_for_Fig3.dat", "GM_fine_for_Fig3.dat", "params_for_coarse_Fig3.dat", and "params_for_fine_Fig3.dat" to generate Figure 3(a) and 3(b).

## Experimental Data
The experimetal data is stored in "Experimental_data", which have two directories storing the experimental data corresponding to Fig.4 and Fig.5 of the paper. 
