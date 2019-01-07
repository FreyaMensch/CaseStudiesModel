%% function to calculate n_h-hourly mean values from hourly measurement data
% in the present case it calculates 6-hourly means in the following ranges:
% 3:00-9:00, 9:00-15:00, 15:00-21:00, 21:00-3:00

function parMean = daily_mean(parameter,n_h) 
%n_b = 24/n_h;               % nr of elements a day is subdivided into
n_b=1;
n_d = floor(length(parameter(n_b:end))/n_h);   % nr of complete days starting at 3:00 on the first day (takes the length of the parameter vector and rounds it down to a full day)
%n_a = n_d*4;                                % multiplied by 4 to provide 4 values per day = length of vector
parMean=zeros(n_d,1);
for i =1:n_d
    parMean(i,1) = mean(parameter(((i-1)*n_h+1):(i*n_h)));
end
end


%% To be changed in the main code:
% meterological input data - daily mean
% n_h = 6;                    % nr of hours that shall be averaged over
% n_b = 24/n_h;               % nr of elements a day is subdivided into
% tm_n=floor(length(t_m));                % time in days [d] but as hourly values 
% tm_n=(floor(length(t_m(n_b:end))/n_h));  % 6-hourly time as fractions of complete days [d]
% t_m=(1:tm_n).*(60*60*24);            % time in seconds [s]
