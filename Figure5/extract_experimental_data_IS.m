%extract data from image to calculate paprams
clear all; close all; clc;
image_IS = "segmented_image_IS.png";
I = imread(image_IS);
greyI = rgb2gray(I);
fiber_idx = find(greyI ==29);
[row_fiber,col_fiber] = ind2sub(size(greyI),fiber_idx);
fiber_row_candidates = unique(row_fiber);
fiber_row_final = fiber_row_candidates(1);
fiber_row_bot = fiber_row_candidates(2);
fiber_thickness = fiber_row_bot- fiber_row_final;
curve_idx = find(greyI >29);

[row_curve,col_curve] = ind2sub(size(greyI),curve_idx);
meaningful_idx = (row_curve < fiber_row_final);
meaningful_row_curve=row_curve(meaningful_idx);
meaningful_col_curve=col_curve(meaningful_idx);
min_fiber_thick = 79;
cropped_image= greyI(min_fiber_thick:fiber_row_final,1:end-1);

%extract desired curve by getting only the minimum of each col
[row_bead,col_bead] = ind2sub(size(cropped_image),find(cropped_image));

%plot image above

extracted_curve = zeros(1,size(cropped_image,2));
for i = 1:size(cropped_image,2)
    ith_col = find(col_bead==i);
    ith_row = min(row_bead(ith_col));
    extracted_curve(i) = ith_row;
end
extracted_curve_IS = size(cropped_image,1)*ones(1,size(extracted_curve,2))-extracted_curve;


%find a scaling in experimental data
dim_fiber_thickness = 0.2;
scale = dim_fiber_thickness/fiber_thickness;
experiment_x_IS = (1:size(extracted_curve_IS,2))* scale;
experiment_u_IS = extracted_curve_IS* scale;
%plot(experiment_x_IS,experiment_u_IS)


% make x,u dimensionless
xscale_IS = 0.0008804057511*1000;
yscale_IS =  0.0003096261559*1000;
nd_x_IS = experiment_x_IS/xscale_IS;
nd_u_IS = experiment_u_IS/yscale_IS;
nd_u_IS_smooth = smooth(nd_u_IS);
for I = 1:2
    nd_u_IS_smooth  = smooth(nd_u_IS_smooth );
end

% find L by getting the first and the second peak
first_peak = (4.80837+4.27833)/2;
first_peak_idx = 120;
second_peak = (44.2978+43.9949)/2;
second_peak_idx = 1166;
second_idx = 0;


% figure;
%  plot(nd_x_IS,nd_u_IS,'k-');
%  hold on;
%  plot(nd_x_IS(first_peak_idx),nd_u_IS(first_peak_idx),'r*');
%  hold on;
%  plot(nd_x_IS(second_peak_idx+second_idx),nd_u_IS(second_peak_idx),'b*');
L_candidate = nd_x_IS(second_peak_idx+second_idx)-nd_x_IS(first_peak_idx);
L = L_candidate;
param_alpha = 3.096261559;
dx = nd_x_IS(2)-nd_x_IS(1);
%calculate numerical_mass
%sanity_check
xx = nd_x_IS(first_peak_idx: second_peak_idx+second_idx)-nd_x_IS(first_peak_idx)*ones(1,length(first_peak_idx: second_peak_idx));
u = nd_u_IS_smooth(first_peak_idx: second_peak_idx+second_idx);
u_experiment = nd_u_IS(first_peak_idx:second_peak_idx+second_idx);
half_alpha_u_sq = param_alpha*0.5*u.*u; 
 q1 = trapz(xx,u);
 q2 = trapz(xx,half_alpha_u_sq);
q_expected = q1+q2;
mass2 = 0;
for i = [first_peak_idx: second_peak_idx+second_idx];
    mass2 = mass2+0.5*dx*(nd_u_IS(i)+nd_u_IS(i+1));
    mass2 = mass2+0.5*dx*(param_alpha*0.5)*(nd_u_IS(i)*nd_u_IS(i)+nd_u_IS(i+1)*nd_u_IS(i+1));
end

Qm = q_expected;
syms x;
eqn = L * (x+ param_alpha * x^2/2)== Qm;
hNs = double(solve(eqn,x));
hN = hNs(2);


%extract data for simulation

%need to make data periodic
u = round(u,4);
idx = find(u==u(1));%1040 is the periodicity thing
idx_for_data = 1:1040;
data = zeros(length(idx_for_data),2);
data(:,1) = xx(idx_for_data);
data(:,2) = u(idx_for_data);

%plot(data(:,1),data(:,2));

experiment = zeros(length(idx_for_data),2);
experiment(:,1) = data(:,1);
experiment(:,2) = u_experiment(idx_for_data);
%plot(experiment(:,1),experiment(:,2))

filename = ['IS_experiment_data.mat'];
save(filename);
