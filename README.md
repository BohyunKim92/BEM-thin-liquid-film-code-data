# BEM-thin-liquid-film-code-data
Accompanying Matlab and c code for "A positivity-preserving numerical method for a thin liquid film on a vertical cylindrical fiber".

B. Kim, H. Ji, A. L. Bertozzi, A. Sadeghpour, and Y. S. Ju. *A positivity-
preserving numerical method for a thin liquid film on a vertical cylindrical fiber*
Under Review in Journal of Computational Physics (2023)

## Algorithms
### Bounded Entropy Method(BEM)
Code for Bounded Entropy Method (Positivity-preserving)

### Generic Method(GM)
Code for Generic Method

## Figures
### Figure 2 data and generating code
Data and code for generating Figure 2 are stored in "Figure2". 
- Running "plot_GM_BEM_Fig2_dimensionless.m" will use "Figure2(a)_GM_coarse_data.dat", "Figure2(b)_BEM_coarse_data.dat", and "parameters_used_for_Fig2.dat" to generate Figure 2(a) and 2(b).

### Figure 3 data and generating code
Data and code for generating Figure 3 are stored in "Figure3". 
- Running "plot_GM_BEM_Fig3_dimensionless.m" will use "BEM_coarse_for_Fig3.dat", "GM_coarse_for_Fig3.dat", "GM_fine_for_Fig3.dat", "params_for_coarse_Fig3.dat", and "params_for_fine_Fig3.dat" to generate Figure 3(a) and 3(b).

### Figure 4 data and generating code
Data and code for generating Figure 4 are stored in "Figure4". 
- Running "create_segmented_experimental_image_RP.m" will preprocess "Experimental_data_RP.png" to produce a segmented image.
- "compare_with_experiment_RP.m" is used to plot experimental data against simulation data using "segmented_RPimage_with_fiber.png" (experimental data), "Figure4_BEM_data.dat", "Figure4_GM_data.dat", and "parameters_used_for_Fig4.dat".

### Figure 5 data and generating code
Data and code for generating Figure 5 are stored in "Figure5". 
- "create_segmented_experimental_image_IS.m" will preprocess "Experimental_data_IS.jpg" to produce a segmented image.
- "extract_experimental_data_IS.m" will then use "segmented_image_IS.png" to extract dat experimental data of IS scheme and save it as
- "compare_with_experiment_IS.m" is used to plot experimental data against simulation data using "IS_experiment_data.dat" (experimental data), "Figure5_BEM_data.dat", "Figure5_GM_data.dat", and "parameters_used_for_Fig5.dat".

### Figure 6 data and generating code
Data and code for generating Figure 6 are stored in "Figure6". 
- "adaptive_time_no_sing.m" is used to plot change of time \delta_t data against simulation time data using "my_pde_adaptive_time_no_sing.dat", "cpu_time_pps_adaptive_time_no_sing.dat" , and "params_adaptive_time_no_sing.dat".

### Figure 7 data and generating code
Data and code for generating Figure 7 and Figure 8 are stored in "Figure7_and_8". 
- "adaptive_singular.m" is used to plot change of time \delta_t data against simulation time data using "my_pde_adaptive_time_no_sing.dat", "cpu_time_pps_adaptive_time_no_sing.dat" , and "params_adaptive_time_no_sing.dat".


