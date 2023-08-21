%adaptive time figure
%singular behavior
%close all; clear all; clc_sing_behavior

myfile = 'my_pde_adaptive_time_singular.dat';
mytime = 'cpu_time_pps_adaptive_time_singular.dat';
params_loc = 'params_adaptive_time_singular.dat'; 
my_pde = importdata(myfile);
params = importdata(params_loc);
my_time = importdata(mytime);

N = params(1,5);
L = params(1,6);
dx = L/(N+1);
it = length(my_pde(:,1))/(N+1);

Position = [100 100 700 500];

% checking time
tts = my_time(:,1);
dts = my_time(:,2);
tts = tts-dts;
my_cpu_time = my_time(:,3); %cpu time for each newton iteration
mytot_calc_time = my_time(:,4); % total calculation time upto here

%% Figure 5.3 global
figure;
indx = tts <1;
%indx = tts<100000;
tts_less_1 = tts(indx);
dts_less_1 = dts(indx);
plot(tts_less_1,dts_less_1,'k');
set(gcf,'position',Position)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XTick = [0 0.25 0.5 0.75 1.0];
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca,'FontSize',30)
xlabel('t')
ylabel('$\Delta t$','Interpreter','latex')


%% Figure 5.3 zoomed
figure;
indx = tts<0.1;
tts_less_25 = tts(indx);
dts_less_25 = dts(indx);
plot(tts_less_25,dts_less_25,'k');
set(gcf,'position',Position)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XTick = [0 0.025 0.05 0.075 0.1];
ax.YTick = [0 10^-3 2.5*10^-3 5.0*10^-3 7.5*10^-3 0.01];
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca,'FontSize',30)
xlabel('t')
ylabel('$\Delta t$','Interpreter','latex')

%plotting when tt = 0.045288 which is around the time of singularity.
%picking interesting points 
% j = [1,10,24,30,34,50]


figure;
tt_sing = [];


it = length(my_pde(:,1))/(N+1);
j = [20,22,25,30, 70];
for i = 1:length(j)
    i = j(i);
    u0_beg_indx = (i-1)*(N+1)+1;
    u0_end_indx =(i-1)*(N+1)+N+1;

    my_x = my_pde(u0_beg_indx:u0_end_indx,1);
    my_u = my_pde(u0_beg_indx:u0_end_indx,2);
    tt_sing = [tt_sing; my_pde(u0_beg_indx,3)];
    if i == 20 
        plot(my_x,my_u,'k:'); hold on;
    elseif i == 22
        plot(my_x,my_u,'k-.'); hold on;
    elseif i == 25
            plot(my_x,my_u,'k--'); hold on;
    elseif i ==30
        indx = 1:4:N+1; % plot every 4 point.
         scatter(my_x(indx),my_u(indx), 'k', 'diamond', 'filled');
         
    elseif i ==70
        plot(my_x,my_u,'k-'); hold on;
    
    end
end
legend(strcat('t=',num2str(tt_sing)))

set(gcf,'position',Position)
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
%ax.XTick = [0 0.025 0.05 0.075 0.1];
%ax.YTick = [0 10^-3 2.5*10^-3 5.0*10^-3 7.5*10^-3 0.01]
set(findall(gca, 'Type', 'Line'),'LineWidth',3)
set(gca,'FontSize',30)
xlabel('x')
ylabel('h')
ylim([-0.1,max(2.5)])
%yticks([10^-3])
%xticklabels(num2str([10^-3]));




%title("tt vs dt");
% figure;
% plot(tts,my_cpu_time);
% title("tt vs cpu time per newton");
% figure;
% plot(tts,mytot_calc_time);
% title("tt vs total cpu time");