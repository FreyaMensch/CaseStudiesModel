%% Plots for temperature model
%% Set interesting time and space variables
% temporal
t_day_s=(10:60:330)*(24*60);      % day of the year to plot [s]  % the indexing is in minutes because the values are every 60sec!
t_day_e=t_day_s+(24*60);          % end of that day
t_day_m=t_day_s+(12*60);          % middle of that day

% spatial
%x_day = ((1:2:reach_nr*2).*(L_reach)/2)./(dx);      % takes the middle of each reach
x_day=[100 121 141]; %.*((L_reach)/2)./(dx); 
x_interval_m=7500;                  % choose spatial interval [m]
x_interval=x_interval_m/dx;         % translates the interval from meters to n-th cell in the vector
[length_matrix, width_matrix]= size(Tw);
nx_pos=floor(width_matrix/x_interval);  % nr of positions calculated (in distance x_interval to eachother)

%% average Tw to smooth curve
n_h = 28*24*60;                     % nr of mins that shall be averaged over
T_avg14d = ones(floor(length_matrix/n_h),floor(width_matrix/x_interval));

figure
for i=1:nx_pos
    j=i*x_interval;
    Tw_avg14d=daily_mean(Tw(:,j),n_h)-273.15;
    T_avg14d(:,i) = Tw_avg14d;
    t_avg=daily_mean(t,n_h);    
    
    plot(t_avg/24/60/60,T_avg14d(:,i))
    hold on
    Legend1{i}=strcat('x=',num2str((x_interval_m*i)/1000),'km'); % along profile');
end
xlabel('Time t [d]')
ylabel('Temperature T [\circC]')
grid on
xlim([0 330])
legend(Legend1); %,'Location','southwest');
title(legend(Legend1),'Position')
%saveas(gcf,[pwd '/figures/Tw_anualTimeSeries_14dayAvg_source.png']);


%% dayily plots: in one plot: different days at fixed location
% plot over one day every 3 months% at three different locations:
% generic: locations and day in the year
figure
for j=1:reach_nr
subplot(1,3,j)
%f(j)=figure
for i=1:length(t_day_s)
    plot((t(1:24*60+1))/60/60,Tw((t_day_s(i):t_day_e(i)),x_day(j))-273.15)
    hold on
    Legend2{i}=strcat(num2str(t_day_s(i)*dt/24/3600),'th day');
end
title(sprintf('position%3.0f km',(x(x_day(j))+dx/2)/1000))
    xticks(0:4:24)
    grid on
    ylim([0 25])
    xlim([0 24])
end
legend(Legend2,'Location','northwest');
set(get(subplot(1,3,1),'YLabel'),'String','Water Temperature T_w [\circC]');
set(get(subplot(1,3,2),'XLabel'),'String','Time t [h]');

%saveas(gcf,[pwd '/figures/Tw_dailyTimeSeries_6days_3locations_source.png']);


%% daily plots: in one plot: different locations at same day 
figure
for i=1:length(t_day_s)
    subplot(1,length(t_day_s),i)
for j=1:reach_nr    
    plot((t(1:24*60+1))/60/60,Tw((t_day_s(i):t_day_e(i)),x_day(j))-273.15)
    hold on
    Legend3{j}=strcat(num2str((x(x_day(j))+dx/2)/1000),'km');
end
    hold on
    plot((t(1:24*60+1))/60/60,T_a_m((t_day_s(i):t_day_e(i)))-273.15,'--')
    %hold off
    title(sprintf('day %3.0f',t_day_s(i)/24/60))
    xticks(0:8:24)
    grid on
    ylim([-10 30])
    xlim([0 24])
end
Legend3{reach_nr+1}=strcat('Air Temperature');
legend(Legend3(1:reach_nr+1),'Location','northwest');
set(get(subplot(1,length(t_day_s),1),'YLabel'),'String','Water Temperature T_w [\circC]');
set(get(subplot(1,length(t_day_s),length(t_day_s)/2),'XLabel'),'String','Time t [h]');

%saveas(gcf,[pwd '/figures/Tw_dailyTimeSeries_3daysAt3locaPlusTa_source.png']);



%% from Vero
figure(6)
for i=1:length(t_day_s)
plot (x/1000,Tw(t_day_m(i),:)-273.15,'-')
hold on
%plot (x/1000,T_a_m(t_day_m(i))-273.15,'--')

Legend4{i}=sprintf('%6.0fth day at %3.0fC air temp',t_day_m(i)*dt/24/3600,T_a_m(t_day_m(i))-273.15);
end
legend(Legend4,'Location','northwest');
xlabel('x [Km]');
ylabel('T [{\circ}C]');
%title('Space variation of Temperature at fixed time')
%saveas(gcf,[pwd '/figures/Tw_spatialVariation_3days10to310_withTa_source.png']);


%% plot meteorological data
%close all
a = 1;             % spacing

% plots solar radiation & air temperature over time in 1plot with 2 y-axis
figure
yyaxis left
plot((t(1:a:end)/24/60/60), H_G_m(1:a:end))
ylabel('solar radiation [W/m^2]')
yyaxis right
plot((t(1:a:end)/24/60/60),T_a_m(1:a:end)-273.15)
ylabel('air temperature [\circC]')
xlabel('time [d]')
%saveas(gcf,[pwd '/figures/input_Ta_HG_year.png']);

 % goal: spatial variation over at one specific day