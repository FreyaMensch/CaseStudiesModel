% function to calculate daily mean values from hourly measurement data

function parMean = daily_mean(parameter) 
n = floor(length(parameter)/24);      
parMean=zeros(n,1);
for i =1:n
    parMean(i,1) = mean(parameter(((i-1)*24+1):(i*24)));
end
end