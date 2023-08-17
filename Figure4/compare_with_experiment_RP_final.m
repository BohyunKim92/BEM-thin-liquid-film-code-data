% scale numerical data to an experiment comparable data
clear all; close all; clc;
myfile = 'sim_data/my_pde_experiment_RP_figure.dat';
hjfile = 'sim_data/out_pde_experiment_RP_figure.dat';
params_loc = 'sim_data/params_experiment_RP_figure.dat';
[parentdir,~,~]=fileparts(pwd);
my_pde = importdata(fullfile(parentdir,myfile));
generic_pde = importdata(fullfile(parentdir,hjfile));
params = importdata(fullfile(parentdir,params_loc));
param_alpha = params(1,1);
param_eta = params(1,2);
param_A = params(1,3);
param_lambda = params(1,4);
N = params(1,5);
L =params(1,6);
dx = params(1,7);
dt0 = params(1,8);
err_tol = params(1,9);
N2 = params(1,10);
N_gen = N;
N_pps = N;
it = length(generic_pde(:,1))/N+1;
k = 12;
figure_size = [18	269	1535	137];


tt = zeros(length(k),1);

for j = 1:length(k)
    i = k(j);
    begindx = (i-1)*(N_gen+1)+1;
    endindx = (i-1)*(N_gen+1)+N_gen+1;
    mybegindx = (i-1)*(N_pps+1)+1;
    myendindx = (i-1)*(N_pps+1)+N_pps+1;
    diff = N_gen/N_pps;
    u_indx = begindx:diff:endindx;
    u_x = generic_pde(u_indx,1);
    u = generic_pde(u_indx,2);
    my_x = my_pde(mybegindx:myendindx,1);
    my_u = my_pde(mybegindx:myendindx,2);
    
end
%<<<<<<< HEAD
gen_tt = generic_pde(u_indx(1),4)
my_tt = my_pde(mybegindx,4)

xscale = 0.001090602322*1000;
yscale =  0.0005885566642*1000;
tscale = 0.001618264896;
gen_tt_dim = gen_tt *tscale
my_tt_dim = my_tt* tscale
%=======

xscale = 0.001090602322*1000;
yscale =  0.0005885566642*1000;
%>>>>>>> 93bebdef967ed4062f11bb875eca151c8e03d397

u_x_scale = u_x*xscale;
u_scale = u*yscale;
my_x_scale =  my_x*xscale; 
my_u_scale = my_u*yscale;
%plot(u_x_scale,u_scale); hold on; plot(my_x_scale,my_u_scale,'--');

%modify data so that it will have L = 30

imax = 7;
Lend = L*xscale*imax;
dx_scale = dx*xscale;
%newN = ceil(newL/dx_scale);
new_x = 0:dx_scale:Lend; 
new_u_pps = [];
new_u_gen = [];
for i = 1:imax
    new_u_pps = [my_u_scale;new_u_pps];
    new_u_gen = [u_scale;new_u_gen];
end
%new_x = new_x(1:newN);
%new_u_gen = new_u_gen(1:newN);
%new_u_pps = new_u_pps(1:newN);

%f = figure;
%f.Position = figure_size;
%plot(new_x,new_u_gen); hold on; 
%plot(new_x,new_u_pps,'--');


%ax = gca;
%ax.FontSize = 25; ax.TickDir = 'out';
%xlabel('x(mm)');
%set(findall(ax, 'Type','Line'),'LineWidth',2); 
%legend("t= "+num2str(tt));
%plot(u_x,u); 
%ylim([-0.8 0])
%xlim([0 30])
%exportgraphics(f,'figures/our_experiment-4.png','Resolution',300)
%saveas(f,'f
%daspect([1 1 1])
% I = imread('figures/our_experiment-4.png');
% crop_position = [75 68 130 112];
% [J,rect] = imcrop(I,crop_position);
% imshow(I);
%imwrite(f,'figures/our_experiment-4.png');

%% plotting experiment image 
image_file = 'segmented_image_with_fiber.png';
I = imread(image_file);
greyI = rgb2gray(I);
fiber_idx = find(greyI ==15);
[row_fiber,col_fiber] = ind2sub(size(greyI),fiber_idx);
fiber_row_candidates=unique(row_fiber);
%choose 134 as the fiber_row
curve_idx = find(greyI ==255);
%imshow(greyI ==255)
[row_curve,col_curve] = ind2sub(size(greyI),curve_idx);
meaningful_idx = (row_curve > 134);
meaningful_row_curve=row_curve(meaningful_idx);
meaningful_col_curve=col_curve(meaningful_idx);
%min_fiber_thick = min(meaningful_row_curve);
max_fiber_thick = 180;%max(meaningful_row_curve);
cropped_image= greyI(134:max_fiber_thick,:);
%imshow(cropped_image)
%extract desired curve
[row_bead,col_bead] = ind2sub(size(cropped_image),find(cropped_image==255));
[C,i_sorted,ic] = unique(col_bead); 
sorted_row_bead = row_bead(i_sorted);
extracted_curve = (sorted_row_bead-ones(size(sorted_row_bead)));

%convert extracted curve in pixel to mm
experiment_u = (extracted_curve * 0.2645833333);
experiment_u = experiment_u';
experiment_x = (1:size(extracted_curve,1))* 0.2645833333;

%%
max_extract = max(extracted_curve)*0.2645833333;
xmax = size(extracted_curve,1)* 0.2645833333;
experiment_u = (extracted_curve * 0.2645833333)*(max(new_u_pps)/max_extract); % scale it same as the pps
experiment_x = (1:size(extracted_curve,1))* 0.2645833333*(max(new_x)/xmax);

%figure_indx = 0:dx_scale:L*xscale*4;
f= figure;
f.Position = [18         355        1239         182]
plot(experiment_x,experiment_u,'k');
hold on;
plot(new_x(1:end-1),new_u_gen,'r--'); hold on; 
hold on;
plot(new_x(1:end-1),new_u_pps,'b:')

%<<<<<<< HEAD
%compare l2 err
%l2_GM = norm()
%l2_BEM = 
%=======
%>>>>>>> 93bebdef967ed4062f11bb875eca151c8e03d397
ax = gca;
ax.FontSize = 20; ax.TickDir = 'out';
xlabel('x(mm)');ylabel('h(mm)')
set(findall(ax, 'Type','Line'),'LineWidth',3);
%legend("Experiment", "BEM",'Location','northeast')
%lgd = legend("Experiment", "BEM",'Location','northeast')
lgd = legend("Experiment","GM", "BEM",'Location','northeast')
%lgd.FontSize = 12;
%ylim([-0.04,1.2]);
xlim([0 22])
%daspect([1 1 1])
%set(gca,'ytick',[min(new_u_pps),max(new_u_pps)])
%set(gca,'yticklabel',{round(min(new_u_pps),3),round(max(new_u_pps),3)})