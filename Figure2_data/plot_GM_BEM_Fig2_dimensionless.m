% compare simulation of two numerical methods on dimensionless scale
% Create figure 2 (a) and 2(b)
% Compare performance of BEM and GM on a coarse grid at physically
% meaningful setting

clear all; clc; close all;
%% loading the data from sim_data folder
hjfile = 'GM_coarse_for_Fig2.dat';
myfile = 'BEM_coarse_for_Fig2.dat';

params_loc = 'params_for_Fig2.dat'; 
generic_pde = importdata(hjfile);
my_pde = importdata(myfile);
params = importdata(params_loc);
%[parentdir,~,~]=fileparts(pwd);
% generic_pde = importdata(fullfile(parentdir,hjfile));
% my_pde = importdata(fullfile(parentdir,myfile));
% params = importdata(fullfile(parentdir,params_loc));

%% extracting parameters, grid size
N = params(1,5);
N_gen = N;
N_pps = N;
L0 = params(1,6);
dx_gen = L0/N_gen;

U0_pps_max = 0; %keep track of maximum value of GM
U0_gen_max = 0; %keep track of maximum value of BEM
U0_pps_min = 4; %keep track of maximum value of GM
U0_gen_min = 4; %keep track of maximum value of BEM

%% initializing vectors for generating figures
k = 1; 
U2_gen = zeros(k*(N+1),1);
U2_pps = zeros(k*(N+1),1);
X2 = zeros(k*(N+1),1);
t2_gen = zeros(k*(N+1),1);
iter = 0;
tts_gen = zeros(5,1); %store each time of GM
tts_pps = zeros(5,1); % store each time of BEM
f = figure(1); % Figure 2a
f1 = figure(2); % Figure 2b
Position = [100 100 1600 1000]; % adjusting figure size& position

%% extracting data and generating 3D plots
for j = 1:2:9
    iter = iter+1;
    u0_beg_indx = (j-1)*(N+1)+1;
    u0_end_indx =(j-1)*(N+1)+N+1;
    tt_gen= generic_pde(u0_beg_indx,4);
    tt_pps= my_pde(u0_beg_indx,4);
    tts_gen(iter) = tt_gen;
    tts_pps(iter) = tt_pps;

    X0_gen = generic_pde(u0_beg_indx:u0_end_indx,1); % x grid
    U0_gen = generic_pde(u0_beg_indx:u0_end_indx,2); % GM data
    U0_pps = my_pde(u0_beg_indx:u0_end_indx,2); % BEM data
    tt_gen = tt_gen*ones(size(U0_gen));
    tt_pps = tt_pps*ones(size(U0_pps));
    
    % make the plot smooth by not plotting all values 
    indx = 1:2:length(X0_gen);
    indx2 = find(U0_gen(indx) <0);
    tt_coarse_gen = tt_gen(indx);
    tt_coarse_pps = tt_pps(indx);
    X0_coarse = X0_gen(indx); % used for both BEM and GM
    U0_coarse_gen = U0_gen(indx);
    U0_coarse_pps = U0_pps(indx);

    % updating min and max value
    if U0_gen_max < max(U0_gen)
        U0_gen_max = max(U0_gen);
    end

    if U0_pps_max < max(U0_pps)
        U0_pps_max = max(U0_pps);
    end
    
    if U0_gen_min > min(U0_gen)
        U0_gen_min = min(U0_gen);
    end

    if U0_pps_min > min(U0_pps)
        U0_pps_min = min(U0_pps);
    end

    figure(1);
    plot3(tt_coarse_gen,X0_coarse,U0_coarse_gen,'k-'); hold on; % plot GM
    figure(1);
    plot3(tt_coarse_gen(indx2),X0_coarse(indx2),U0_coarse_gen(indx2),'bs','MarkerSize',40); hold on; %plot negative values and indicate the location with blue square 
    figure(2);
    plot3(tt_coarse_pps,X0_coarse,U0_coarse_pps,'k-'); hold on; % plot BEM
end

angle = 190+68;

%% adjust figures for GM Fig 2(a)
figure(1);
ax = gca;
ax.XGrid = 'on';
ax.ZGrid = 'on';
ax.YMinorGrid = 'on';
set(gcf,'position',Position)
set(gca,'ydir','reverse')
set(gca,'LineWidth',2);
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca,'FontSize',40)
xticks(tts_gen)
xticklabels(num2str(round(tts_gen,3)));
zticks([floor(U0_gen_min)+1:ceil(U0_gen_max)])
xlabel('t')
ylabel('x')
zlabel('h')
xlim([floor(tts_gen(1)) tts_gen(end)])
ylim([0 L0])
zlim([(U0_gen_min-0.1),ceil(U0_gen_max)])
view(angle,25)
pbaspect([12 15 3])

%% adjust figures for BEM Fig 2(b)
figure(2);
ax = gca;
ax.XGrid = 'on';
ax.ZGrid = 'on';
ax.YMinorGrid = 'on';
set(gcf,'position',Position)
set(gca,'ydir','reverse')
set(gca,'LineWidth',2);
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca,'FontSize',40)
xticks(tts_pps)
xticklabels(num2str(round(tts_pps,3)));
zticks([floor(U0_pps_min) floor(U0_pps_min)+1:ceil(U0_pps_max)])
xlabel('t')
ylabel('x ')
zlabel('h')
xlim([floor(tts_pps(1)) tts_pps(end)])
ylim([0 L0])
zlim([floor(U0_pps_min), ceil(U0_pps_max)])
view(angle,25)
pbaspect([12 15 3])


