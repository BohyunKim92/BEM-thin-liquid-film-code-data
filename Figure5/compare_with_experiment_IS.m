% scale numerical data to an experiment comparable data
%clear all; close all; clc;
load("IS_experiment_data.mat");
myfile = 'Figure5_BEM_data.dat';
hjfile = 'Figure5_GM_data.dat';
params_loc = 'parameters_used_for_Fig5.dat';

[parentdir,~,~]=fileparts(pwd);
my_pde = importdata(myfile);
generic_pde = importdata(hjfile);
params = importdata(params_loc);
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
k_out = 73;
k_my = 80;

xscale_IS = 0.0008804057511*1000;
yscale_IS =  0.0003096261559*1000;
tscale = 0.002581268242;

experiment_x_scale = experiment(:,1)* xscale_IS;
experiment_u_scale = experiment(:,2)* yscale_IS;
% extract data from generic scheme
i_out = k_out;
begindx = (i_out-1)*(N_gen+1)+1;
endindx = (i_out-1)*(N_gen+1)+N_gen+1;
diff = N_gen/N_pps;
u_indx = begindx:diff:endindx;
u_x = generic_pde(u_indx,1);
u = generic_pde(u_indx,2);
tt_u = generic_pde(u_indx(1),4);
tt_u_dimensional = tt_u * tscale;

u_x_scale = u_x*xscale_IS;
u_scale = u*yscale_IS ;


f= figure;
f.Position = [18         355        1239         182];
%do circular shif to make the data come to the middle
experiment_u_scale = circshift(experiment_u_scale,length(experiment_u_scale)/2);
plot(experiment_x_scale,experiment_u_scale,'k');
hold on;
u_scale = circshift(u_scale,length(u_scale)/2);

plot(u_x_scale ,u_scale, 'r--'); hold on; 



%extract data from BEM scheme
i_my = k_my;
mybegindx = (i_my-1)*(N_pps+1)+1;
myendindx = (i_my-1)*(N_pps+1)+N_pps+1;
my_x = my_pde(mybegindx:myendindx,1);
my_u = my_pde(mybegindx:myendindx,2);
my_tt = my_pde(mybegindx,4);
my_tt_dimensional = my_tt*tscale;

my_x_scale =  my_x*xscale_IS; 
my_u_scale = my_u*yscale_IS;
hold on;
my_u_scale = circshift(my_u_scale,length(my_u_scale)/2);
plot(my_x_scale,my_u_scale,'b:')


ax = gca;
ax.FontSize = 20; ax.TickDir = 'out';
xlabel('x(mm)');ylabel('h(mm)')
set(findall(ax, 'Type','Line'),'LineWidth',3);
lgd = legend("Experiment","GM", "BEM",'Location','northeast');
xlim([experiment_x_scale(1,1), experiment_x_scale(end,1)])