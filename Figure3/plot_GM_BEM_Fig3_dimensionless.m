% compare simulation of two numerical methods on dimensionless scale
% Create figure 3 (a) and 3(b)
% Compare performance of BEM and GM on a coarse grid at physically
% meaningful setting
% Snap shot of figure 2

clear all; clc; close all;


%% loading the data from sim_data folder
mycoarsefile = 'BEM_coarse_for_Fig3.dat';
hjcoarsefile = 'GM_coarse_for_Fig3.dat';
hjfinefile = 'GM_fine_for_Fig3.dat';
params_loc = 'params_for_coarse_Fig3.dat'; 
params_loc2 = 'params_for_fine_Fig3.dat'; 
my_pde_coarse = importdata(mycoarsefile);
generic_pde_coarse = importdata(hjcoarsefile);
generic_pde_fine= importdata(hjfinefile);

%% extracting parameters, grid size
params = importdata(params_loc);
params2 = importdata(params_loc2);
param_alpha = params(1,1);
param_eta = params(1,2);
param_A = params(1,3);
param_lambda = params(1,4);
N_coarse =params(1,5);
N_fine = params2(1,5);
L =params(1,6);
dt0 = params(1,8);
xx = params(2:end,1);
u = params(2:end,2);


%% setting up the figure and time
figure_size1 = [100 100 1200 700]; % figure size
figure_size2 = [100 100 1000 600]; % figure size

time = 145; % capture the profile at this time step
k_gen_c = [time];
k_gen_f = [time];
k_pps_c = [time];

f1 = figure;
f2 = figure;
f1.Position = figure_size1;
f2.Position = figure_size2;
dx = [L/N_coarse; L/N_coarse; L/N_fine]; 


%% plotting figures
for j = 1:length(k_gen_c)
   i_gen_c = k_gen_c(j);
   i_gen_f = k_gen_f(j);
   i_pps_c = k_pps_c(j);

   begindx_coarse = (i_gen_c-1)*(N_coarse+1)+1;
   endindx_coarse = (i_gen_c-1)*(N_coarse+1)+N_coarse+1;
   mybegindx_coarse = (i_pps_c-1)*(N_coarse+1)+1;
   myendindx_coarse = (i_pps_c-1)*(N_coarse+1)+N_coarse+1;
   begindx_fine = (i_gen_f-1)*(N_fine+1)+1;
   endindx_fine = (i_gen_f-1)*(N_fine+1)+N_fine+1;

   u_x_coarse = generic_pde_coarse(begindx_coarse:endindx_coarse,1);
   u_x_fine = generic_pde_fine(begindx_fine:endindx_fine,1);
   u_coarse = generic_pde_coarse(begindx_coarse:endindx_coarse,2);
   u_fine = generic_pde_fine(begindx_fine:endindx_fine,2);
   my_x = my_pde_coarse(mybegindx_coarse:myendindx_coarse,1);
   my_u = my_pde_coarse(mybegindx_coarse:myendindx_coarse,2);
%    
%    u_x_coarse = generic_pde_coarse(begindx_coarse:endindx_coarse,1)*xscale;
%    u_x_fine = generic_pde_fine(begindx_fine:endindx_fine,1)*xscale;
%    u_coarse = generic_pde_coarse(begindx_coarse:endindx_coarse,2)*yscale;
%    u_fine = generic_pde_fine(begindx_fine:endindx_fine,2)*yscale;
%    my_x = my_pde_coarse(mybegindx_coarse:myendindx_coarse,1)*xscale;
%    my_u = my_pde_coarse(mybegindx_coarse:myendindx_coarse,2)*yscale;
  
   tt_coarse = generic_pde_coarse(begindx_coarse,4)
   tt_fine= generic_pde_fine(begindx_fine,4)
   mytt_coarse = my_pde_coarse(mybegindx_coarse,4)
   [M,I] = min(u_coarse);
   % only plot a few for smoothness
   smoothindx_c = 1:2:N_coarse+1;
   smoothindx_f = 1:4:N_fine+1;

   figure(1)
   plot(u_x_coarse(smoothindx_c),u_coarse(smoothindx_c),'r-.'); hold on;plot(my_x(smoothindx_c),my_u(smoothindx_c),'b:'); hold on; plot(u_x_fine(smoothindx_f), u_fine(smoothindx_f),'k-'); hold on;
   plot(u_x_coarse(I), M, 'bs','MarkerSize',40)
   figure(2)
   plot(u_x_coarse(smoothindx_c),u_coarse(smoothindx_c),'r-.'); hold on;plot(my_x(smoothindx_c),my_u(smoothindx_c),'b:'); hold on; plot(u_x_fine(smoothindx_f), u_fine(smoothindx_f),'k-'); hold on;
   plot(u_x_coarse(I), M, 'bs','MarkerSize',40); hold on;
   plot(u_x_coarse(smoothindx_c),zeros(length(u_x_coarse(smoothindx_c))),'k--')  
end

%compare l2 err
indx = 1:2:length(u_x_fine);
figure;
plot(u_x_coarse, u_fine(indx)-u_coarse);hold on;
plot(u_x_coarse, u_fine(indx)-my_u);
err_GM = norm(u_fine(indx)-u_coarse,2)/(L);
err_BEM = norm(u_fine(indx)-my_u,2)/(L);
legend('GM err', 'BEM err')

figure(1)
ax = gca;
ax.FontSize = 25; ax.TickDir = 'out';
xlabel('x');ylabel('h')
set(findall(ax, 'Type','Line'),'LineWidth',4);
legend("Coarse GM simulation","Coarse BEM simulation", "Fine GM simulation",'Location','northwest')
ylim([min(u_coarse),ceil(max(u_coarse))]);
xlim([0 L])
set(gca,'xtick',[0,5,10,15,20,24])
daspect([2.0 1 1])

figure(2)
ax = gca;
ax.FontSize = 25; ax.TickDir = 'out';
xlabel('x');ylabel('h')
set(findall(ax, 'Type','Line'),'LineWidth',4);
legend("Coarse GM simulation","Coarse BEM simulation", "Fine GM simulation",'Location','northwest')
ylim([min(u_coarse),0.07]);
xlim([19.0 L])
set(gca,'ytick',[min(u_coarse),0.0,0.07])
%set(gca,'yticklabel',{num2str(round(min(u_coarse),4)),'0','0.05'})
%daspect([1 1 1])